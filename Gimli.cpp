#include "Gimli.h"
#include <iostream>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <iomanip>
#include <string>
#include <fstream>
#include <algorithm>
#include "gurobi_c++.h"
using namespace std;

bool compare(Column c0, Column c1) {
	if (c0.y > c1.y) {
		return true;
	}
	else if (c0.y == c1.y) {
		if (c0.z > c1.z) {
			return true;
		}
	}
	return false;
}

bool compareX(Column c0, Column c1) {
	if (c0.x > c1.x) {
		return true;
	}
	return false;
}

Gimli::Gimli() {
	//c=0x9e377900
	UINT32 value = 0;
	for (int i = 24; i > 0; i--) {
		if (i % 4 == 0) {
			value = 0x9e377900 ^ i;
			for (int j = 0; j < 32; j++) {
				constant[6 - i / 24][j] = (value >> j) & 0x1;
			}
		}
	}
	threeVarLinearZero = new int* [4];
	for (int i = 0; i < 4; i++) {
		threeVarLinearZero[i] = new int[5];
	}

	andConstraint = new int* [9];
	for (int i = 0; i < 9; i++) {
		andConstraint[i] = new int[7];
	}

	fourVarLinearZero = new int* [8];
	for (int i = 0; i < 8; i++) {
		fourVarLinearZero[i] = new int[6];
	}

	orConstraint = new int* [9];
	for (int i = 0; i < 9; i++) {
		orConstraint[i] = new int[7];
	}

	andValueConstraint = new int* [12];
	for (int i = 0; i < 12; i++) {
		andValueConstraint[i] = new int[7];
	}

	orValueConstraint = new int* [12];
	for (int i = 0; i < 12; i++) {
		orValueConstraint[i] = new int[7];
	}

	//initializeMinimalModel(".\\constraint\\threeLinearZero.txt", threeVarLinearZero);
	//initializeMinimalModel(".\\constraint\\and.txt", andConstraint);
	//initializeMinimalModel(".\\constraint\\fourLinearZero.txt", fourVarLinearZero);
	//initializeMinimalModel(".\\constraint\\or.txt", orConstraint);
	//initializeMinimalModel(".\\constraint\\andValue.txt", andValueConstraint);
	//initializeMinimalModel(".\\constraint\\orValue.txt", orValueConstraint);

	initializeMinimalModel("threeLinearZero.txt", threeVarLinearZero);
	initializeMinimalModel("and.txt", andConstraint);
	initializeMinimalModel("fourLinearZero.txt", fourVarLinearZero);
	initializeMinimalModel("or.txt", orConstraint);
	initializeMinimalModel("andValue.txt", andValueConstraint);
	initializeMinimalModel("orValue.txt", orValueConstraint);

	srand(time(NULL));
	previous = 0;
	previousINT32 = 0;
	solNum = 0;
}

Gimli::~Gimli() {
	for (int i = 0; i < 4; i++) {
		delete[]threeVarLinearZero[i];
	}
	for (int i = 0; i < 9; i++) {
		delete[]andConstraint[i];
		delete[]orConstraint[i];
	}
	for (int i = 0; i < 8; i++) {
		delete[]fourVarLinearZero[i];
	}
	for (int i = 0; i < 12; i++) {
		delete[]andValueConstraint[i];
		delete[]orValueConstraint[i];
	}
	delete[]threeVarLinearZero;
	delete[]andConstraint;
	delete[]orConstraint;
	delete[]fourVarLinearZero;
	delete[]andValueConstraint;
	delete[]orValueConstraint;
}

UINT64 Gimli::getRand64() {
	UINT64 sum = rand() % EXP[16];
	sum = (sum << 16) | (rand() % EXP[16]);
	sum = (sum << 16) | (rand() % EXP[16]);
	sum = (sum << 16) | (rand() % EXP[16]);
	sum = sum + previous;
	previous = sum;
	return sum;
}

UINT32 Gimli::getRand32() {
	UINT32 sum = rand() % EXP[16];
	sum = (sum << 16) | (rand() % EXP[16]);
	sum = sum + previousINT32;
	previousINT32 = sum;
	return sum;
}

void Gimli::initializeMinimalModel(string filename,int **matrix) {
	int row, col;
	ifstream inFile(filename);
	inFile >> row >> col;
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			inFile >> matrix[i][j];
		}
	}
	inFile.close();
}

UINT32 Gimli::LL(UINT32 number, int n) {
	if (n == 0) {
		return number;
	}
	return (number << n) | (number >> (32 - n));
}

UINT32 Gimli::RR(UINT32 number, int n) {
	if (n == 0) {
		return number;
	}
	return (number >> n) | (number << (32 - n));
}

void Gimli::SPBox(UINT32 *input, UINT32 *output) {
	UINT32 x = LL(input[0], 24);
	UINT32 y = LL(input[1], 9);
	UINT32 z = input[2];

	output[2] = x ^ (z << 1) ^ ((y & z) << 2);
	output[1] = y ^ x ^ ((x | z) << 1);
	output[0] = z ^ y ^ ((x & y) << 3);
}

void Gimli::inverseSPBox(UINT32 input[]) {
	UINT32 tx = input[0], ty = input[1], tz = input[2];
	bool xBit[32], yBit[32], zBit[32];
	for (int i = 0; i < 32; i++) {
		xBit[i] = (tx >> i) & 0x1;
		yBit[i] = (ty >> i) & 0x1;
		zBit[i] = (tz >> i) & 0x1;
	}
	//bit level computation
	bool reXBit[32], reYBit[32], reZBit[32];
	//Firstly, compute reXBit[0],reYBit[0],reZBit[0]
	reXBit[0] = zBit[0];
	reYBit[0] = yBit[0] ^ reXBit[0];
	reZBit[0] = xBit[0] ^ reYBit[0];
	//Secondly, compute reXBit[1],reYBit[1],reZBit[1]
	reXBit[1] = zBit[1] ^ reZBit[0];
	reYBit[1] = yBit[1] ^ reXBit[1] ^ (reXBit[0] | reZBit[0]);
	reZBit[1] = xBit[1] ^ reYBit[1];
	//Thirdly, compute reXBit[2],reYBit[2],reZBit[2]
	reXBit[2] = zBit[2] ^ reZBit[1] ^ (reYBit[0] & reZBit[0]);
	reYBit[2] = yBit[2] ^ reXBit[2] ^ (reXBit[1] | reZBit[1]);
	reZBit[2] = xBit[2] ^ reYBit[2];
	//Recursively recover
	for (int i = 3; i < 32; i++) {
		reXBit[i] = zBit[i] ^ reZBit[i - 1] ^ (reYBit[i - 2] & reZBit[i - 2]);
		reYBit[i] = yBit[i] ^ reXBit[i] ^ (reXBit[i - 1] | reZBit[i - 1]);
		reZBit[i] = xBit[i] ^ reYBit[i] ^ (reXBit[i - 3] & reYBit[i - 3]);
	}
	UINT32 x = 0;
	UINT32 y = 0;
	UINT32 z = 0;
	for (int i = 0; i < 32; i++) {
		if (reZBit[i]) {
			z |= EXP[i];
		}
		if (reYBit[i]) {
			y |= EXP[(i - 9 + 32) % 32];
		}
		if (reXBit[i]) {
			x |= EXP[(i - 24 + 32) % 32];
		}
	}
	input[0] = x;
	input[1] = y;
	input[2] = z;
}

void Gimli::modelXUpdate(GRBModel& model, vector<GRBVar>& id, vector<GRBVar>& od, vector<GRBVar>& odn, vector<GRBVar>& v) {
	//id: input difference
	//od: output difference
	//odn: output difference of nonlinear part

	// x = x<<<24
	// y = y<<<9
	// x' = z + y + (x & y)<<3

	//First condier x'[0,1,2]
	//Model the linear relation (three variables)
	//in[0] + in[1] - out >= 0
	//in[0] - in[1] + out >= 0
	//-in[0] + in[1] + out >= 0
	//-in[0] - in[1] - out + 2 >= 0

	//x'[i] = z[i] + y[i-9]
	for (int i = 0; i < 3; i++) {
		loadConstraintThreeLinearZero(id[ZMap(i)], id[YMap(MOD(i-9))], od[i], model);
	}

	//consider the remaining bits x'[i] = z[i] + y[i-9] + x[i-27]y[i-12]
	for (int i = 3; i < 32; i++) {
		loadAndConstraint(v[MOD(i - 27)], v[YMap(MOD(i - 12))], id[MOD(i - 27)], id[YMap(MOD(i - 12))], odn[i], model);
		//load the relation between (iz[i],iy[i-9],odn[i],od[i])//three variable linear
		loadConstraintFourLinearZero(id[ZMap(i)], id[YMap(MOD(i - 9))], odn[i], od[i], model);
	}
}

void Gimli::modelYUpdate(GRBModel& model, vector<GRBVar>& id, vector<GRBVar>& od, vector<GRBVar>& odn, vector<GRBVar>& v) {
	//id: input difference
	//od: output difference
	//odn: output difference of nonlinear part

	// x = x<<<24
	// y = y<<<9
	// y' = y + x + (x | z)<<1

	//Firstly, consider y'[0]
	//y'[0] = y[0-9] + x[0-24]
	loadConstraintThreeLinearZero(id[YMap(MOD(0 - 9))], id[MOD(0 - 24)], od[YMap(0)], model);

	//Secondly, consider y'[i] (i>=1)
	//y'[i] = y[i-9] + x[i-24] + (x[i-25] | z[i-1])
	for (int i = 1; i < 32; i++) {
		loadORConstraint(v[MOD(i - 25)], v[ZMap(MOD(i - 1))], id[MOD(i - 25)], id[ZMap(MOD(i - 1))], odn[YMap(i)], model);
		//load linear relation of difference
		loadConstraintFourLinearZero(id[YMap(MOD(i - 9))], id[MOD(i - 24)], odn[YMap(i)], od[YMap(i)], model);
	}
}

void Gimli::modelZUpdate(GRBModel& model, vector<GRBVar>& id, vector<GRBVar>& od, vector<GRBVar>& odn, vector<GRBVar>& v) {
	//id: input difference
	//od: output difference
	//odn: output difference of nonlinear part

	// x = x<<<24
	// y = y<<<9
	// z' = x + z<<1 + (y & z)<<2

	//z'[0] = x[0-24]
	model.addConstr(id[MOD(0 - 24)] == od[ZMap(0)]);

	//z'[1] = x[1-24] + z[0]
	loadConstraintThreeLinearZero(id[MOD(1 - 24)], id[ZMap(0)], od[ZMap(1)], model);

	//z'[i] = x[i-24] + z[i-1] + y[i-11]&z[i-2]
	for (int i = 2; i < 32; i++) {
		loadAndConstraint(v[YMap(MOD(i - 11))], v[ZMap(MOD(i - 2))], id[YMap(MOD(i - 11))], id[ZMap(MOD(i - 2))], odn[ZMap(i)], model);
		loadConstraintFourLinearZero(id[MOD(i - 24)], id[ZMap(MOD(i - 1))], odn[ZMap(i)], od[ZMap(i)], model);
	}
}

void Gimli::modelXUpdate(GRBModel& model, vector<GRBVar>& idX, vector<GRBVar>& idY, vector<GRBVar>& idZ, vector<GRBVar>& odX, vector<GRBVar>& odnX, vector<GRBVar>& vX, vector<GRBVar>& vY) {
	// x = x<<<24
	// y = y<<<9
	// x' = z + y + (x & y)<<3
	//x'[i] = z[i] + y[i-9]
	for (int i = 0; i < 3; i++) {
		loadConstraintThreeLinearZero(idZ[i], idY[MOD(i - 9)], odX[i], model);
	}

	//consider the remaining bits x'[i] = z[i] + y[i-9] + x[i-27]y[i-12]
	for (int i = 3; i < 32; i++) {
		loadAndConstraint(vX[MOD(i - 27)], vY[MOD(i - 12)], idX[MOD(i - 27)], idY[MOD(i - 12)], odnX[i], model);
		//load the relation between (iz[i],iy[i-9],odn[i],od[i])//three variable linear
		loadConstraintFourLinearZero(idZ[i], idY[MOD(i - 9)], odnX[i], odX[i], model);
	}
}

void Gimli::modelYUpdate(GRBModel& model, vector<GRBVar>& idX, vector<GRBVar>& idY, vector<GRBVar>& idZ, vector<GRBVar>& odY, vector<GRBVar>& odnY, vector<GRBVar>& vX, vector<GRBVar>& vZ) {
	// x = x<<<24
	// y = y<<<9
	// y' = y + x + (x | z)<<1

	//Firstly, consider y'[0]
	//y'[0] = y[0-9] + x[0-24]
	loadConstraintThreeLinearZero(idY[MOD(0 - 9)], idX[MOD(0 - 24)], odY[0], model);

	//Secondly, consider y'[i] (i>=1)
	//y'[i] = y[i-9] + x[i-24] + (x[i-25] | z[i-1])
	for (int i = 1; i < 32; i++) {
		loadORConstraint(vX[MOD(i - 25)], vZ[MOD(i - 1)], idX[MOD(i - 25)], idZ[MOD(i - 1)], odnY[i], model);
		//load linear relation of difference
		loadConstraintFourLinearZero(idY[MOD(i - 9)], idX[MOD(i - 24)], odnY[i], odY[i], model);
	}
}

void Gimli::modelZUpdate(GRBModel& model, vector<GRBVar>& idX, vector<GRBVar>& idY, vector<GRBVar>& idZ, vector<GRBVar>& odZ, vector<GRBVar>& odnZ, vector<GRBVar>& vY, vector<GRBVar>& vZ) {
	// x = x<<<24
	// y = y<<<9
	// z' = x + z<<1 + (y & z)<<2

	//z'[0] = x[0-24]
	model.addConstr(idX[MOD(0 - 24)] == odZ[0]);

	//z'[1] = x[1-24] + z[0]
	loadConstraintThreeLinearZero(idX[MOD(1 - 24)], idZ[0], odZ[1], model);

	//z'[i] = x[i-24] + z[i-1] + y[i-11]&z[i-2]
	for (int i = 2; i < 32; i++) {
		loadAndConstraint(vY[MOD(i - 11)], vZ[MOD(i - 2)], idY[MOD(i - 11)], idZ[MOD(i - 2)], odnZ[i], model);
		loadConstraintFourLinearZero(idX[MOD(i - 24)], idZ[MOD(i - 1)], odnZ[i], odZ[i], model);
	}
}

void Gimli::loadConstraintThreeLinearZero(GRBVar& a, GRBVar& b, GRBVar& c, GRBModel& model) {
	for (int i = 0; i < 4; i++) {
		model.addConstr(threeVarLinearZero[i][0] * a + threeVarLinearZero[i][1] * b + threeVarLinearZero[i][2] * c + threeVarLinearZero[i][3] >= threeVarLinearZero[i][4]);
	}
}

void Gimli::loadConstraintFourLinearZero(GRBVar& a, GRBVar& b, GRBVar& c, GRBVar& d, GRBModel& model) {
	for (int i = 0; i < 8; i++) {
		model.addConstr(fourVarLinearZero[i][0] * a + fourVarLinearZero[i][1] * b + fourVarLinearZero[i][2] * c + fourVarLinearZero[i][3]*d + fourVarLinearZero[i][4]>= fourVarLinearZero[i][5]);
	}
}

void Gimli::loadAndConstraint(GRBVar& a, GRBVar& b, GRBVar& c, GRBVar& d, GRBVar& e, GRBModel& model) {
	for (int i = 0; i < 9; i++) {
		model.addConstr(andConstraint[i][0] * a + andConstraint[i][1] * b + andConstraint[i][2] * c + andConstraint[i][3] * d + andConstraint[i][4] * e + andConstraint[i][5] >= andConstraint[i][6]);
	}
}

void Gimli::loadORConstraint(GRBVar& a, GRBVar& b, GRBVar& c, GRBVar& d, GRBVar& e, GRBModel& model) {
	for (int i = 0; i < 9; i++) {
		model.addConstr(orConstraint[i][0] * a + orConstraint[i][1] * b + orConstraint[i][2] * c + orConstraint[i][3] * d + orConstraint[i][4] * e + orConstraint[i][5] >= orConstraint[i][6]);
	}
}

void Gimli::loadConstraintOnTheOutputDifference(GRBModel& model, vector<GRBVar>& di, UINT32 output[]) {
	//output[0] = x -> di[0-31]
	//output[1] = y -> di[32-63]
	//output[2] = z -> di[64-95]
	int value = 0;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 32; j++) {
			value = (output[i] >> j) & 0x1;
			model.addConstr(di[i * 32 + j] == value);
		}
	}
}

void Gimli::loadConstraintOnTheOutputSingleDifference(GRBModel& model, vector<GRBVar>& di, UINT32 output) {
	int value = 0;
	for (int j = 0; j < 32; j++) {
		value = (output >> j) & 0x1;
		model.addConstr(di[j] == value);
	}
}

void Gimli::modelXValueUpdate(GRBModel& model, vector<GRBVar>& iv, vector<GRBVar>& ov) {
	// x = x<<<24
	// y = y<<<9
	// x' = z + y + (x & y)<<3

	//First consider x'[0,1,2]
	//x'[i] = z[i] + y[i-9]
	for (int i = 0; i < 3; i++) {
		loadConstraintThreeLinearZero(iv[ZMap(i)], iv[YMap(MOD(i - 9))], ov[i], model);
	}

	//Next consider the nonlinear part
	//x'[i] = z[i] + y[i-9] + x[i-27]&y[i-12]
	for (int i = 3; i < 32; i++) {
		loadAndValueConstraint(iv[ZMap(i)], iv[YMap(MOD(i - 9))], iv[MOD(i-27)], iv[YMap(MOD(i-12))], ov[i],model);
	}
}

void Gimli::modelYValueUpdate(GRBModel& model, vector<GRBVar>& iv, vector<GRBVar>& ov) {
	// x = x<<<24
	// y = y<<<9
	// y' = y + x + (x | z)<<1

	//Firstly, consider y'[0]
	//y'[0] = y[0-9] + x[0-24]
	loadConstraintThreeLinearZero(iv[YMap(MOD(0 - 9))], iv[MOD(0 - 24)], ov[YMap(0)], model);

	//Secondly, consider y'[i] (i>=1)
	//y'[i] = y[i-9] + x[i-24] + (x[i-25] | z[i-1])
	for (int i = 1; i < 32; i++) {
		loadORValueConstraint(iv[YMap(MOD(i - 9))], iv[MOD(i - 24)], iv[MOD(i - 25)], iv[ZMap(MOD(i - 1))], ov[YMap(i)], model);
	}
}

void Gimli::modelZValueUpdate(GRBModel& model, vector<GRBVar>& iv, vector<GRBVar>& ov) {
	// x = x<<<24
	// y = y<<<9
	// z' = x + z<<1 + (y & z)<<2

	//z'[0] = x[0-24]
	model.addConstr(iv[MOD(0 - 24)] == ov[ZMap(0)]);

	//z'[1] = x[1-24] + z[0]
	loadConstraintThreeLinearZero(iv[MOD(1 - 24)], iv[ZMap(0)], ov[ZMap(1)], model);

	//z'[i] = x[i-24] + z[i-1] + y[i-11]&z[i-2]
	for (int i = 2; i < 32; i++) {
		loadAndValueConstraint(iv[MOD(i - 24)], iv[ZMap(MOD(i - 1))], iv[YMap(MOD(i - 11))], iv[ZMap(MOD(i - 2))], ov[ZMap(i)], model);
	}
}

void Gimli::modelXValueUpdate(GRBModel& model, vector<GRBVar>& iX, vector<GRBVar>& iY, vector<GRBVar>& iZ, vector<GRBVar>& oX) {
	// x = x<<<24
	// y = y<<<9
	// x' = z + y + (x & y)<<3

	//First consider x'[0,1,2]
	//x'[i] = z[i] + y[i-9]
	for (int i = 0; i < 3; i++) {
		//loadConstraintThreeLinearZero(iZ[i], iY[MOD(i - 9)], oX[i], model);
		/*
		1 1 -1 0 0
		1 -1 1 0 0
		-1 1 1 0 0
		-1 -1 -1 2 0
		*/
		model.addConstr(iZ[i] + iY[MOD(i - 9)] - oX[i] >= 0);
		model.addConstr(iZ[i] - iY[MOD(i - 9)] + oX[i] >= 0);
		model.addConstr(-iZ[i] + iY[MOD(i - 9)] + oX[i] >= 0);
		model.addConstr(-iZ[i] - iY[MOD(i - 9)] - oX[i] + 2 >= 0);
	}

	//Next consider the nonlinear part
	//x'[i] = z[i] + y[i-9] + x[i-27]&y[i-12]
	for (int i = 3; i < 32; i++) {
		/*
		-1 1 0 1 1 0 0
		1 -1 0 1 1 0 0
		1 1 1 0 -1 0 0==
		1 1 0 1 -1 0 0==
		1 -1 1 0 1 0 0==
		-1 1 1 0 1 0 0==
		1 1 -1 -1 1 1 0==
		-1 -1 1 0 -1 2 0==
		1 -1 -1 -1 -1 3 0==
		-1 -1 -1 -1 1 3 0==
		-1 -1 0 1 -1 2 0==
		-1 1 -1 -1 -1 3 0
		*/
		//loadAndValueConstraint(iZ[i], iY[MOD(i - 9)], iX[MOD(i - 27)], iY[MOD(i - 12)], oX[i], model);
		model.addConstr(-iZ[i] + iY[MOD(i - 9)] + iY[MOD(i - 12)] + oX[i] >= 0);
		model.addConstr(iZ[i] - iY[MOD(i - 9)] + iY[MOD(i - 12)] + oX[i] >= 0);
		model.addConstr(iZ[i] + iY[MOD(i - 9)] + iX[MOD(i - 27)] - oX[i] >= 0);
		model.addConstr(iZ[i] + iY[MOD(i - 9)] + iY[MOD(i - 12)] - oX[i] >= 0);
		model.addConstr(iZ[i] - iY[MOD(i - 9)] + iX[MOD(i - 27)] + oX[i] >= 0);
		model.addConstr(-iZ[i] + iY[MOD(i - 9)] + iX[MOD(i - 27)] + oX[i] >= 0);
		model.addConstr(iZ[i] + iY[MOD(i - 9)] - iX[MOD(i - 27)] - iY[MOD(i - 12)] + oX[i] + 1 >= 0);
		model.addConstr(-iZ[i] - iY[MOD(i - 9)] + iX[MOD(i - 27)] - oX[i] + 2 >= 0);
		model.addConstr(iZ[i] - iY[MOD(i - 9)] - iX[MOD(i - 27)] - iY[MOD(i - 12)] - oX[i] + 3 >= 0);
		model.addConstr(-iZ[i] - iY[MOD(i - 9)] - iX[MOD(i - 27)] - iY[MOD(i - 12)] + oX[i] + 3 >= 0);
		model.addConstr(-iZ[i] - iY[MOD(i - 9)] + iY[MOD(i - 12)] - oX[i] + 2 >= 0);
		model.addConstr(-iZ[i] + iY[MOD(i - 9)] - iX[MOD(i - 27)] - iY[MOD(i - 12)] - oX[i] + 3 >= 0);
	}
}

void Gimli::modelYValueUpdate(GRBModel& model, vector<GRBVar>& iX, vector<GRBVar>& iY, vector<GRBVar>& iZ, vector<GRBVar>& oY) {
	// x = x<<<24
	// y = y<<<9
	// y' = y + x + (x | z)<<1

	//Firstly, consider y'[0]
	//y'[0] = y[0-9] + x[0-24]
	loadConstraintThreeLinearZero(iY[MOD(0 - 9)], iX[MOD(0 - 24)], oY[0], model);

	//Secondly, consider y'[i] (i>=1)
	//y'[i] = y[i-9] + x[i-24] + (x[i-25] | z[i-1])
	for (int i = 1; i < 32; i++) {
		loadORValueConstraint(iY[MOD(i - 9)], iX[MOD(i - 24)], iX[MOD(i - 25)], iZ[MOD(i - 1)], oY[i], model);
	}
}

void Gimli::modelZValueUpdate(GRBModel& model, vector<GRBVar>& iX, vector<GRBVar>& iY, vector<GRBVar>& iZ, vector<GRBVar>& oZ) {
	// x = x<<<24
	// y = y<<<9
	// z' = x + z<<1 + (y & z)<<2

	//z'[0] = x[0-24]
	model.addConstr(iX[MOD(0 - 24)] == oZ[0]);

	//z'[1] = x[1-24] + z[0]
	//loadConstraintThreeLinearZero(iX[MOD(1 - 24)], iZ[0], oZ[1], model);
	/*
		1 1 -1 0 0
		1 -1 1 0 0
		-1 1 1 0 0
		-1 -1 -1 2 0
		*/
	model.addConstr(iX[MOD(1 - 24)]+iZ[0]-oZ[1] >= 0);
	model.addConstr(iX[MOD(1 - 24)]-iZ[0]+oZ[1] >= 0);
	model.addConstr(-iX[MOD(1 - 24)]+iZ[0]+oZ[1] >= 0);
	model.addConstr(-iX[MOD(1 - 24)]-iZ[0]-oZ[1] + 2 >= 0);

	//z'[i] = x[i-24] + z[i-1] + y[i-11]&z[i-2]
	for (int i = 2; i < 32; i++) {
		//loadAndValueConstraint(iX[MOD(i - 24)], iZ[MOD(i - 1)], iY[MOD(i - 11)], iZ[MOD(i - 2)], oZ[i], model);
		/*
		-1 1 0 1 1 0 0
		1 -1 0 1 1 0 0==
		1 1 1 0 -1 0 0==
		1 1 0 1 -1 0 0==
		1 -1 1 0 1 0 0==
		-1 1 1 0 1 0 0==
		1 1 -1 -1 1 1 0==
		-1 -1 1 0 -1 2 0==
		1 -1 -1 -1 -1 3 0==
		-1 -1 -1 -1 1 3 0==
		-1 -1 0 1 -1 2 0==
		-1 1 -1 -1 -1 3 0
		*/
		model.addConstr(-iX[MOD(i - 24)]+ iZ[MOD(i - 1)]+iZ[MOD(i - 2)]+ oZ[i] >= 0);
		model.addConstr(iX[MOD(i - 24)]- iZ[MOD(i - 1)]+ iZ[MOD(i - 2)]+ oZ[i] >= 0);
		model.addConstr(iX[MOD(i - 24)]+ iZ[MOD(i - 1)]+ iY[MOD(i - 11)]- oZ[i] >= 0);
		model.addConstr(iX[MOD(i - 24)]+ iZ[MOD(i - 1)]+ iZ[MOD(i - 2)]- oZ[i] >= 0);
		model.addConstr(iX[MOD(i - 24)]- iZ[MOD(i - 1)]+ iY[MOD(i - 11)]+ oZ[i] >= 0);
		model.addConstr(-iX[MOD(i - 24)]+ iZ[MOD(i - 1)]+ iY[MOD(i - 11)]+oZ[i] >= 0);
		model.addConstr(iX[MOD(i - 24)]+ iZ[MOD(i - 1)]- iY[MOD(i - 11)]- iZ[MOD(i - 2)]+ oZ[i] + 1 >= 0);
		model.addConstr(-iX[MOD(i - 24)]- iZ[MOD(i - 1)]+ iY[MOD(i - 11)]- oZ[i] + 2 >= 0);
		model.addConstr(iX[MOD(i - 24)]- iZ[MOD(i - 1)]- iY[MOD(i - 11)]- iZ[MOD(i - 2)]- oZ[i] + 3 >= 0);
		model.addConstr(-iX[MOD(i - 24)]- iZ[MOD(i - 1)]- iY[MOD(i - 11)]- iZ[MOD(i - 2)]+ oZ[i] + 3 >= 0);
		model.addConstr(-iX[MOD(i - 24)]- iZ[MOD(i - 1)]+ iZ[MOD(i - 2)]- oZ[i] + 2 >= 0);
		model.addConstr(-iX[MOD(i - 24)]+ iZ[MOD(i - 1)]- iY[MOD(i - 11)]- iZ[MOD(i - 2)]- oZ[i] + 3 >= 0);
	}
}

void Gimli::loadAndValueConstraint(GRBVar& a, GRBVar& b, GRBVar& c, GRBVar& d, GRBVar& e, GRBModel& model) {
	for (int i = 0; i < 12; i++) {
		model.addConstr(andValueConstraint[i][0] * a + andValueConstraint[i][1] * b + andValueConstraint[i][2] * c + andValueConstraint[i][3] * d + andValueConstraint[i][4] * e + andValueConstraint[i][5] >= andValueConstraint[i][6]);
	}
}

void Gimli::loadORValueConstraint(GRBVar& a, GRBVar& b, GRBVar& c, GRBVar& d, GRBVar& e, GRBModel& model) {
	for (int i = 0; i < 12; i++) {
		model.addConstr(orValueConstraint[i][0] * a + orValueConstraint[i][1] * b + orValueConstraint[i][2] * c + orValueConstraint[i][3] * d + orValueConstraint[i][4] * e + orValueConstraint[i][5] >= orValueConstraint[i][6]);
	}
}

void Gimli::wordToBit(UINT32 word, bool* bitArray) {
	for (int i = 0; i < 32; i++) {
		bitArray[i] = (word >> i) & 0x1;
	}
}

void Gimli::extraLinearRelationsFromX(int bitPos, bool* IX, bool* IY, bool* IZ, bool* OX, int pos[], int coef[], bool& value) {
	// x = x<<<24
	// y = y<<<9
	// x' = z + y + (x & y)<<3

	//First consider x'[0,1,2]
	//x'[i] = z[i] + y[i-9]

	bool nd = OX[bitPos] ^ IZ[bitPos] ^ IY[MOD(bitPos - 9)];//The difference of the nonlinear part
	int xPos= MOD(bitPos - 27);
	int yPos = MOD(bitPos - 12);
	pos[0] = MOD(bitPos - 27);
	pos[1] = YMap(MOD(bitPos - 12));
	coef[0] = 0;
	coef[1] = 0;
	value = 0;
	if (bitPos < 3) {
		if (nd) {
			cout << "Contradition X in " << bitPos << endl;
		}
	}
	else {
		if (nd && IX[xPos] && IY[yPos]) {//1 1 1 (ab=c)==> a=b
			coef[0] = 1;
			coef[1] = 1;
			value = 0;
		}
		else if (nd && IX[xPos] && !IY[yPos]) {//1 0 1 (ab=c)==> b=1
			coef[1] = 1;
			value = 1;
		}
		else if (nd && !IX[xPos] && IY[yPos]) {//0 1 1 (ab=c)==> a=1
			coef[0] = 1;
			value = 1;
		}
		else if (nd && !IX[xPos] && !IY[yPos]) {//0 0 1 (ab=c)==> contradiction
			cout << "Contradition X in " << bitPos << endl;
		}
		//
		else if (!nd && IX[xPos] && IY[yPos]) {//1 1 0 (ab=c)==> a=b
			coef[0] = 1;
			coef[1] = 1;
			value = 1;
		}
		else if (!nd && IX[xPos] && !IY[yPos]) {//1 0 0 (ab=c)==> b=0
			coef[1] = 1;
			value = 0;
		}
		else if (!nd && !IX[xPos] && IY[yPos]) {//0 1 0 (ab=c)==> a=0
			coef[0] = 1;
			value = 0;
		}
	}
}

void Gimli::extraLinearRelationsFromY(int bitPos, bool* IX, bool* IY, bool* IZ, bool* OY, int pos[], int coef[], bool& value) {
	// x = x<<<24
	// y = y<<<9
	// y' = y + x + (x | z)<<1

	//Firstly, consider y'[0]
	//y'[0] = y[0-9] + x[0-24]

	//Secondly, consider y'[i] (i>=1)
	//y'[i] = y[i-9] + x[i-24] + (x[i-25] | z[i-1])

	bool nd = OY[bitPos] ^ IY[MOD(bitPos-9)] ^ IX[MOD(bitPos-24)];//The difference of the nonlinear part
	int xPos = MOD(bitPos - 25);
	int zPos = MOD(bitPos - 1);
	pos[0] = MOD(bitPos - 25);
	pos[1] = ZMap(MOD(bitPos - 1));
	coef[0] = 0;
	coef[1] = 0;
	value = 0;
	if (bitPos < 1) {
		if (nd) {
			cout << "Contradition Y in " << bitPos << endl;
		}
	}
	else {
		if (nd && IX[xPos] && IZ[zPos]) {//1 1 1 a|b=c=> a=b
			coef[0] = 1;
			coef[1] = 1;
			value = 0;
		}
		else if (nd && IX[xPos] && !IZ[zPos]) {//1 0 1 a|b=c=> b=0
			coef[1] = 1;
			value = 0;
		}
		else if (nd && !IX[xPos] && IZ[zPos]) {//0 1 1 a|b=c=> a=0
			coef[0] = 1;
			value = 0;
		}
		else if (nd && !IX[xPos] && !IZ[zPos]) {//0 0 1 a|b=c=> a=0
			cout << "Contradition Y in " << bitPos << endl;
		}
		else if (!nd && IX[xPos] && IZ[zPos]) {//1 1 0 a|b=c=> a+b=1
			coef[0] = 1;
			coef[1] = 1;
			value = 1;
		}
		else if (!nd && IX[xPos] && !IZ[zPos]) {//1 0 0 a|b=c=> b=1
			coef[1] = 1;
			value = 1;
		}
		else if (!nd && !IX[xPos] && IZ[zPos]) {//0 1 1 a|b=c=> a=1
			coef[0] = 1;
			value = 1;
		}
	}
}

void Gimli::extraLinearRelationsFromZ(int bitPos, bool* IX, bool* IY, bool* IZ, bool* OZ, int pos[], int coef[], bool& value) {
	// x = x<<<24
	// y = y<<<9
	// z' = x + z<<1 + (y & z)<<2
	if (bitPos == 0) {
		if (OZ[0] != IX[MOD(0 - 24)]) {
			cout << "Contradition Z in " << bitPos << endl;
		}
		return;
	}

	bool nd = OZ[bitPos] ^ IX[MOD(bitPos-24)] ^ IZ[MOD(bitPos - 1)];//The difference of the nonlinear part
	int yPos = MOD(bitPos - 11);
	int zPos = MOD(bitPos - 2);
	pos[0] = YMap(MOD(bitPos - 11));
	pos[1] = ZMap(MOD(bitPos - 2));
	coef[0] = 0;
	coef[1] = 0;
	value = 0;
	if (bitPos < 2) {
		if (nd) {
			cout << "Contradition Z in " << bitPos << endl;
		}
	}
	else {
		if (nd && IY[yPos] && IZ[zPos]) {//1 1 1 (ab=c)==> a=b
			coef[0] = 1;
			coef[1] = 1;
			value = 0;
		}
		else if (nd && IY[yPos] && !IZ[zPos]) {//1 0 1 (ab=c)==> b=1
			coef[1] = 1;
			value = 1;
		}
		else if (nd && !IY[yPos] && IZ[zPos]) {//0 1 1 (ab=c)==> a=1
			coef[0] = 1;
			value = 1;
		}
		else if (nd && !IY[yPos] && !IZ[zPos]) {//0 0 1 (ab=c)==> contradiction
			cout << "Contradition Z in " << bitPos << endl;
		}
		//
		else if (!nd && IY[yPos] && IZ[zPos]) {//1 1 0 (ab=c)==> a=b
			coef[0] = 1;
			coef[1] = 1;
			value = 1;
		}
		else if (!nd && IY[yPos] && !IZ[zPos]) {//1 0 0 (ab=c)==> b=0
			coef[1] = 1;
			value = 0;
		}
		else if (!nd && !IY[yPos] && IZ[zPos]) {//0 1 0 (ab=c)==> a=0
			coef[0] = 1;
			value = 0;
		}
	}
}

void Gimli::outputCondition(int pos[], int coef[],bool value,int columnNum) {
	int p0 = pos[0] % 32;
	int p1 = pos[1] % 32;
	if (coef[0] || coef[1]) {
		if (coef[0]==1 && coef[1]==1) {
			cout << sign[pos[0] / 32]<<"["<<p0 << "] + " << sign[pos[1] / 32]<< "[" << p1 << "] = " << value;
		}
		else if (coef[0] == 1 && coef[1] == 0) {
			cout << sign[pos[0] / 32]<<"[" << p0 << "] = " << value;
		}
		else if (coef[0] == 0 && coef[1] == 1) {
			cout << sign[pos[1] / 32] << "[" << p1 << "] = " << value;
		}
		cout << endl;
	}
}

int Gimli::addRelationToSystem(int pos[], int coef[], bool& value, bool** system, int row, bool *equation, int eqSize) {
	memset(equation, 0, eqSize);
	if (coef[0]) {
		equation[pos[0]] = 1;
	}
	if (coef[1]){
		equation[pos[1]] = 1;
	}
	equation[eqSize - 1] = value;
	return addToEquationSystem(equation, system, row, eqSize);
}

int Gimli::addToEquationSystem(bool equation[], bool** system, int row, int column) {
	for (int i = 0; i < column-1; i++) {
		if (equation[i]) {//use Gauss elimination on equation
			if (system[i][i]) {
				for (int j = i; j < column; j++) {
					equation[j] ^= system[i][j];
				}
			}
			else {
				//Add to the equation system
				for (int j = i; j < column; j++) {
					system[i][j] = equation[j];
				}
				//Apply Gauss elimination for the i-th row
				for (int j = i + 1; j < row; j++) {
					if (system[j][j] && system[i][j]) {//Add row[j]
						for (int k = j; k < column; k++) {
							system[i][k] ^= system[j][k];
						}
					}
				}
				//Apply Gauss elimination for the i-th row (up)
				for (int j = 0; j < i; j++) {
					if (system[j][j] && system[j][i]) {//Add row[j]
						//Apply gauss elimination to the j-th row (add the i-th row)
						for (int k = i; k < column; k++) {
							system[j][k] ^= system[i][k];
						}
					}
				}
				return 1;//independent
			}
		}
	}
	if (equation[column-1] == 1) {
		//cout << "Contradiction" << endl;
		return 0;
	}
	return 2;//dependent
}

void Gimli::outputMatrix(bool** matrix, int row, int column) {
	bool isFirst = true;
	bool isEmpty = true;
	bool isConstant[192];
	for (int i = 0; i < 192; i++) {
		isConstant[i] = false;
	}

	char valueString[32];
	for (int i = 0; i < row; i++) {
		isFirst = true;
		isEmpty = true;
		for (int j = 0; j < column-1; j++) {
			if (matrix[i][j]) {
				if (isFirst) {
					cout <<j/96<< sign[(j % 96) / 32]<<"["<<j%32<<"]";
					isConstant[i] = true;
					isFirst = false;
					isEmpty = false;
				}
				else {
					cout << " + " << j / 96 << sign[(j % 96) / 32] << "[" << j % 32 << "]";
					isConstant[i] = false;
				}
			}
		}
		if (!isEmpty) {
			cout << " = " << matrix[i][column - 1] << endl;
		}
		/*if (isConstant[i]) {
			if (i % 32 == 0) {
				if (matrix[i][192]) {
					valueString[i%32]='1';
				}
				else {
					valueString[i % 32] = '0';
				}
			}
			else {
				if (matrix[i][192]) {
					valueString[i % 32] = '1';
				}
				else {
					valueString[i % 32] = '0';
				}
			}
		}
		else {
			valueString[i % 32] = '-';
		}

		if ((i+1) % 32 == 0) {
			//output value in binary
			cout << i / 96 << sign[(i % 96) / 32] << ": ";
			for (int k = 31; k >= 0; k--) {
				cout <<"& $"<< valueString[k]<<"$";
				if (k % 4 == 0) {
					cout << " ";
				}
			}
			cout << endl;
		}*/
	}

	//output two-bit conditions
	/*for (int i = 0; i < row; i++) {
		isFirst = true;
		isEmpty = true;
		if (isConstant[i]) {
			continue;
		}
		for (int j = 0; j < column - 1; j++) {
			if (matrix[i][j]) {
				if (isFirst) {
					cout << j / 96 << sign[(j % 96) / 32] << "[" << j % 32 << "]";
					isEmpty = false;
					isFirst = false;
				}
				else {
					cout << " + " << j / 96 << sign[(j % 96) / 32] << "[" << j % 32 << "]";
				}
			}
		}
		if (!isEmpty) {
			cout << " = " << matrix[i][column - 1] << endl;
		}
	}*/
}

string Gimli::toBinary(UINT32 num) {
	string str = "";
	unsigned int left = 0;
	int array[32];
	for (int i = 0; i < 32; i++) {
		array[i] = 0;
	}
	int count = 0;
	while (num) {
		left = num % 2;
		array[count++] = left;
		num = num / 2;
	}

	for (int i = 31; i >= 0; i--) {
		cout << array[i];
		if (array[i] == 1) {
			str += '1';
		}
		else {
			str += '0';
		}
		if (i % 4 == 0) {
			cout << " ";
			str += ' ';
		}
	}
	cout << endl;
	return str;
}

void Gimli::clearMatrix(bool** matrix, int row, int column) {
	for (int i = 0; i < row; i++) {
		memset(matrix[i], 0, column);
	}
}

int Gimli::countFreedom(bool** matrix, int row, int column) {
	int conditionNum = 0;
	for (int i = 0; i < row; i++) {
		if (matrix[i][i]) {
			conditionNum++;
		}
	}
	return conditionNum;
}

void Gimli::matrixEqual(bool** sourceMatrix, bool** desMatrix, int row, int column) {
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < column; j++) {
			desMatrix[i][j] = sourceMatrix[i][j];
		}
	}
}

//verify the 6-round SFS collisions
void Gimli::checkRealSFSCollision() {
	UINT32 state[12] = { 0xff792f16 ,0x37feedd1,0xaca93960,0x9a757bef,0xd8080e8,0x88cda05b,
						0xff792f16 ,0x37feedd1,0xaca93960,0x9a757bef,0xd8080e8,0x88cda05b };
	UINT32 newState[12];
	for (int k = 0; k < 12; k++) {
		newState[k] = state[k];
	}
	newState[3] = state[3] ^ 0x7c2c642a;
	newState[9] = state[9] ^ 0x7c2c642a;

	//outputState(state);
	//cout << "after adding differece:" << endl;
	//outputState(newState);

	//start 5-round permutation
	permutation(state, 6);
	permutation(newState, 6);
	//outputState(state);
	//cout << "after adding differece:" << endl;
	//outputState(newState);

	for (int k = 0; k < 12; k++) {
		newState[k] = newState[k]^state[k];
	}

	cout << "Difference in the output S^6:" << endl;
	outputState(newState);
}

void Gimli::SPColumn(UINT32 state[], int size) {
	UINT32 x, y, z;
	for (int i = 0; i < size / 3; i++) {
		x = LL(state[3 * i], 24);
		y = LL(state[3 * i + 1], 9);
		z = state[3 * i + 2];

		state[3 * i + 2] = x ^ (z << 1) ^ ((y & z) << 2);
		state[3 * i + 1] = y ^ x ^ ((x | z) << 1);
		state[3 * i] = z ^ y ^ ((x & y) << 3);
	}
}

void Gimli::inverseSPColumn(UINT32 state[], int size) {
	UINT32 stateTmp[3];
	for (int i = 0; i < size / 3; i++) {
		stateTmp[0] = state[i * 3];
		stateTmp[1] = state[i * 3+1];
		stateTmp[2] = state[i * 3+2];
		inverseSPBox(stateTmp);

		state[i * 3] = stateTmp[0];
		state[i * 3+1] = stateTmp[1];
		state[i * 3+2] = stateTmp[2];
	}
}

void Gimli::bigSwap(UINT32 state[]) {
	swap(state[0], state[6]);
	swap(state[3], state[9]);
}
void Gimli::smallSwap(UINT32 state[]) {
	swap(state[0], state[3]);
	swap(state[6], state[9]);
}

void Gimli::permutation(UINT32 state[],int rounds) {
	for (int i = 24; i > 24-rounds; i--) {
		SPColumn(state,12);
		if (i % 4 == 0) {//small swap
			swap(state[0], state[3]);
			swap(state[6], state[9]);
		}
		else if (i % 2 == 0) {//Big swap
			swap(state[0], state[6]);
			swap(state[3], state[9]);
		}

		if (i % 4 == 0) {
			state[0] = state[0]^0x9e377900 ^ i;
		}
	}
}

void Gimli::permutation(UINT32 state[], int start, int end) {
	for (int i = start; i > end+1; i--) {
		SPColumn(state, 12);
		if (i % 4 == 0) {//small swap
			swap(state[0], state[3]);
			swap(state[6], state[9]);
		}
		else if (i % 2 == 0) {//Big swap
			swap(state[0], state[6]);
			swap(state[3], state[9]);
		}

		if (i % 4 == 0) {
			state[0] = state[0] ^ 0x9e377900 ^ i;
		}
	}
}

void Gimli::outputState(UINT32 state[]) {
	cout << "State: " << endl;
	for (int r = 0; r < 3; r++) {
		for (int c = 0; c < 4; c++) {
			cout << hex << setw(9) << state[r + 3 * c] << " ";
		}
		cout << endl;
	}
	cout << endl;
}

bool Gimli::testCorrectOfDiff(UINT32 inputDiff[][12], int rounds,int start, int end) {
	GRBEnv env = GRBEnv();
	env.set(GRB_IntParam_OutputFlag, 0);
	env.set(GRB_IntParam_Threads, 3);
	GRBModel model = GRBModel(env);
	vector<vector<GRBVar> > ac;//before constant addition
	vector<vector<vector<GRBVar> > > d;//difference
	vector<vector<vector<GRBVar> > > dn;//difference of nonlinear part
	vector<vector<vector<GRBVar> > > v;//value
	d.resize(rounds);
	dn.resize(rounds);
	v.resize(rounds);
	ac.resize(6);//6 in total
	int variableNum = 32;

	for (int i = 0; i < 6; i++) {
		ac[i].resize(variableNum);
		for (int j = 0; j < 32; j++) {
			ac[i][j]= model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		}
	}

	for (int i = 0; i < rounds; i++) {
		d[i].resize(12);
		v[i].resize(12);
		dn[i].resize(12);
	}
	for (int r = 0; r < rounds; r++) {
		for (int b = 0; b < 12; b++) {
			d[r][b].resize(variableNum);
			v[r][b].resize(variableNum);
			dn[r][b].resize(variableNum);
		}
	}

	for (int r = 0; r < rounds; r++) {
		for (int b = 0; b < 12; b++) {
			for (int z = 0; z < variableNum; z++) {
				dn[r][b][z] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
				d[r][b][z] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
				v[r][b][z] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
			}
		}
	}

	//load diff
	for (int i = 0; i < rounds; i++) {
		for (int j = 0; j < 12; j++) {
			loadConstraintOnTheOutputSingleDifference(model, d[i][j], inputDiff[i][j]);
		}
	}

	int r = 0;
	for (int i = start; i > end+1; i--) {
		r = start - i;
		if (i % 4 == 0) {//Small swap
			modelXUpdate(model, d[r][0], d[r][1], d[r][2], d[r + 1][3], dn[r][0], v[r][0], v[r][1]);//small swap
			modelYUpdate(model, d[r][0], d[r][1], d[r][2], d[r + 1][1], dn[r][1], v[r][0], v[r][2]);
			modelZUpdate(model, d[r][0], d[r][1], d[r][2], d[r + 1][2], dn[r][2], v[r][1], v[r][2]);
			
			modelXUpdate(model, d[r][3], d[r][4], d[r][5], d[r + 1][0], dn[r][3], v[r][3], v[r][4]);
			modelYUpdate(model, d[r][3], d[r][4], d[r][5], d[r + 1][4], dn[r][4], v[r][3], v[r][5]);
			modelZUpdate(model, d[r][3], d[r][4], d[r][5], d[r + 1][5], dn[r][5], v[r][4], v[r][5]);
			
			modelXUpdate(model, d[r][6], d[r][7], d[r][8], d[r + 1][9], dn[r][6], v[r][6], v[r][7]);//small swap
			modelYUpdate(model, d[r][6], d[r][7], d[r][8], d[r + 1][7], dn[r][7], v[r][6], v[r][8]);
			modelZUpdate(model, d[r][6], d[r][7], d[r][8], d[r + 1][8], dn[r][8], v[r][7], v[r][8]);
			
			modelXUpdate(model, d[r][9], d[r][10], d[r][11], d[r + 1][6], dn[r][9], v[r][9], v[r][10]);
			modelYUpdate(model, d[r][9], d[r][10], d[r][11], d[r + 1][10], dn[r][10], v[r][9], v[r][11]);
			modelZUpdate(model, d[r][9], d[r][10], d[r][11], d[r + 1][11], dn[r][11], v[r][10], v[r][11]);
		}
		else if (i % 2 == 0) {//Big swap
			modelXUpdate(model, d[r][0], d[r][1], d[r][2], d[r + 1][6], dn[r][0], v[r][0], v[r][1]);//big swap
			modelYUpdate(model, d[r][0], d[r][1], d[r][2], d[r + 1][1], dn[r][1], v[r][0], v[r][2]);
			modelZUpdate(model, d[r][0], d[r][1], d[r][2], d[r + 1][2], dn[r][2], v[r][1], v[r][2]);

			modelXUpdate(model, d[r][3], d[r][4], d[r][5], d[r + 1][9], dn[r][3], v[r][3], v[r][4]);
			modelYUpdate(model, d[r][3], d[r][4], d[r][5], d[r + 1][4], dn[r][4], v[r][3], v[r][5]);
			modelZUpdate(model, d[r][3], d[r][4], d[r][5], d[r + 1][5], dn[r][5], v[r][4], v[r][5]);

			modelXUpdate(model, d[r][6], d[r][7], d[r][8], d[r + 1][0], dn[r][6], v[r][6], v[r][7]);//big swap
			modelYUpdate(model, d[r][6], d[r][7], d[r][8], d[r + 1][7], dn[r][7], v[r][6], v[r][8]);
			modelZUpdate(model, d[r][6], d[r][7], d[r][8], d[r + 1][8], dn[r][8], v[r][7], v[r][8]);

			modelXUpdate(model, d[r][9], d[r][10], d[r][11], d[r + 1][3], dn[r][9], v[r][9], v[r][10]);
			modelYUpdate(model, d[r][9], d[r][10], d[r][11], d[r + 1][10], dn[r][10], v[r][9], v[r][11]);
			modelZUpdate(model, d[r][9], d[r][10], d[r][11], d[r + 1][11], dn[r][11], v[r][10], v[r][11]);
		}
		else {//No swap
			modelXUpdate(model, d[r][0], d[r][1], d[r][2], d[r + 1][0], dn[r][0], v[r][0], v[r][1]);//big swap
			modelYUpdate(model, d[r][0], d[r][1], d[r][2], d[r + 1][1], dn[r][1], v[r][0], v[r][2]);
			modelZUpdate(model, d[r][0], d[r][1], d[r][2], d[r + 1][2], dn[r][2], v[r][1], v[r][2]);

			modelXUpdate(model, d[r][3], d[r][4], d[r][5], d[r + 1][3], dn[r][3], v[r][3], v[r][4]);
			modelYUpdate(model, d[r][3], d[r][4], d[r][5], d[r + 1][4], dn[r][4], v[r][3], v[r][5]);
			modelZUpdate(model, d[r][3], d[r][4], d[r][5], d[r + 1][5], dn[r][5], v[r][4], v[r][5]);

			modelXUpdate(model, d[r][6], d[r][7], d[r][8], d[r + 1][6], dn[r][6], v[r][6], v[r][7]);//big swap
			modelYUpdate(model, d[r][6], d[r][7], d[r][8], d[r + 1][7], dn[r][7], v[r][6], v[r][8]);
			modelZUpdate(model, d[r][6], d[r][7], d[r][8], d[r + 1][8], dn[r][8], v[r][7], v[r][8]);

			modelXUpdate(model, d[r][9], d[r][10], d[r][11], d[r + 1][9], dn[r][9], v[r][9], v[r][10]);
			modelYUpdate(model, d[r][9], d[r][10], d[r][11], d[r + 1][10], dn[r][10], v[r][9], v[r][11]);
			modelZUpdate(model, d[r][9], d[r][10], d[r][11], d[r + 1][11], dn[r][11], v[r][10], v[r][11]);
		}

		//model value transition
		if (i % 4 == 0) {//Small swap (and constant addition)
			modelXValueUpdate(model, v[r][0], v[r][1], v[r][2], v[r + 1][3]);//small_swap
			modelYValueUpdate(model, v[r][0], v[r][1], v[r][2], v[r + 1][1]);
			modelZValueUpdate(model, v[r][0], v[r][1], v[r][2], v[r + 1][2]);

			//before v[r+1][0], there is a constant addition and we need to consider it.
			//ac[6-i/4]-> v[r+1][0]
			//modelXValueUpdate(model, v[r][3], v[r][4], v[r][5], v[r + 1][0]);//small_swap
			modelXValueUpdate(model, v[r][3], v[r][4], v[r][5], ac[6-i/4]);//small_swap
			//ac -> v[r+1][0]
			for (int j = 0; j < 32; j++) {
				if (constant[6 - i / 4][j]==0) {
					model.addConstr(ac[6 - i / 4][j] == v[r + 1][0][j]);
				}
				else {
					model.addConstr(ac[6 - i / 4][j] + v[r + 1][0][j] == 1);
				}
			}
			modelYValueUpdate(model, v[r][3], v[r][4], v[r][5], v[r + 1][4]);
			modelZValueUpdate(model, v[r][3], v[r][4], v[r][5], v[r + 1][5]);

			modelXValueUpdate(model, v[r][6], v[r][7], v[r][8], v[r + 1][9]);//small_swap
			modelYValueUpdate(model, v[r][6], v[r][7], v[r][8], v[r + 1][7]);
			modelZValueUpdate(model, v[r][6], v[r][7], v[r][8], v[r + 1][8]);

			modelXValueUpdate(model, v[r][9], v[r][10], v[r][11], v[r + 1][6]);//small_swap
			modelYValueUpdate(model, v[r][9], v[r][10], v[r][11], v[r + 1][10]);
			modelZValueUpdate(model, v[r][9], v[r][10], v[r][11], v[r + 1][11]);
		}
		else if (i % 2 == 0) {
			modelXValueUpdate(model, v[r][0], v[r][1], v[r][2], v[r + 1][6]);//small_swap
			modelYValueUpdate(model, v[r][0], v[r][1], v[r][2], v[r + 1][1]);
			modelZValueUpdate(model, v[r][0], v[r][1], v[r][2], v[r + 1][2]);

			modelXValueUpdate(model, v[r][3], v[r][4], v[r][5], v[r + 1][9]);//small_swap
			modelYValueUpdate(model, v[r][3], v[r][4], v[r][5], v[r + 1][4]);
			modelZValueUpdate(model, v[r][3], v[r][4], v[r][5], v[r + 1][5]);

			modelXValueUpdate(model, v[r][6], v[r][7], v[r][8], v[r + 1][0]);//small_swap
			modelYValueUpdate(model, v[r][6], v[r][7], v[r][8], v[r + 1][7]);
			modelZValueUpdate(model, v[r][6], v[r][7], v[r][8], v[r + 1][8]);

			modelXValueUpdate(model, v[r][9], v[r][10], v[r][11], v[r + 1][3]);//small_swap
			modelYValueUpdate(model, v[r][9], v[r][10], v[r][11], v[r + 1][10]);
			modelZValueUpdate(model, v[r][9], v[r][10], v[r][11], v[r + 1][11]);
		}
		else {
			modelXValueUpdate(model, v[r][0], v[r][1], v[r][2], v[r + 1][0]);
			modelYValueUpdate(model, v[r][0], v[r][1], v[r][2], v[r + 1][1]);
			modelZValueUpdate(model, v[r][0], v[r][1], v[r][2], v[r + 1][2]);

			modelXValueUpdate(model, v[r][3], v[r][4], v[r][5], v[r + 1][3]);
			modelYValueUpdate(model, v[r][3], v[r][4], v[r][5], v[r + 1][4]);
			modelZValueUpdate(model, v[r][3], v[r][4], v[r][5], v[r + 1][5]);

			modelXValueUpdate(model, v[r][6], v[r][7], v[r][8], v[r + 1][6]);//(There is a contradiction)
			modelYValueUpdate(model, v[r][6], v[r][7], v[r][8], v[r + 1][7]);
			modelZValueUpdate(model, v[r][6], v[r][7], v[r][8], v[r + 1][8]);

			modelXValueUpdate(model, v[r][9], v[r][10], v[r][11], v[r + 1][9]);
			modelYValueUpdate(model, v[r][9], v[r][10], v[r][11], v[r + 1][10]);
			modelZValueUpdate(model, v[r][9], v[r][10], v[r][11], v[r + 1][11]);
		}
	}
	model.optimize();

	if (model.get(GRB_IntAttr_Status) == 3) {
		cout << "The model is infeasible" << endl;
		return false;
	}

	UINT32 value[12];
	UINT32 state[12];
	cout << "Diff trail:" << endl;
	for (int r = 0; r < rounds; r++) {
		for (int b = 0; b < 12; b++) {
			value[b] = 0;
			for (int i = 0; i < 32; i++) {
				if (d[r][b][i].get(GRB_DoubleAttr_X) == 1) {
					value[b] = value[b] | EXP[i];
				}
			}
		}
		//output value
		for (int i = 0; i < 3; i++) {
			for (int col = 0; col < 4; col++) {
				if (i == 0) {
					//cout << hex << setw(9) << LL(value[i + col * 3],24) << " ";
					cout << hex << setw(9) << value[i + col * 3] << " ";
				}
				else if (i == 1) {
					//cout << hex << setw(9) << LL(value[i + col * 3], 9) << " ";
					cout << hex << setw(9) << value[i + col * 3] << " ";
				}
				else {
					//cout << hex << setw(9) << value[i + col * 3] << " ";
					cout << hex << setw(9) << value[i + col * 3] << " ";
				}
			}
			cout << endl;
		}
		cout << endl;
	}

	cout << "Message satisfying the differential conditions:" << endl;
	for (int r = 0; r < rounds; r++) {
		for (int b = 0; b < 12; b++) {
			value[b] = 0;
			for (int i = 0; i < 32; i++) {
				if (v[r][b][i].get(GRB_DoubleAttr_X) == 1) {
					value[b] = value[b] | EXP[i];
				}
			}
		}
		if (r == 0) {
			for (int i = 0; i < 12; i++) {
				state[i] = value[i];
			}
		}
		//output value
		for (int i = 0; i < 3; i++) {
			for (int col = 0; col < 4; col++) {
				cout << hex << setw(9) << value[i + col * 3] << " ";
			}
			cout << endl;
		}
		cout << endl;
	}

	UINT32 newState[12];
	for (int i = 0; i < 12; i++) {
		newState[i] = state[i] ^ inputDiff[0][i];
	}
	permutation(state, start, end);
	permutation(newState, start, end);

	for (int i = 0; i < 12; i++) {
		state[i] = state[i] ^ newState[i];
	}

	cout << "Final output difference (verification of the results returned by the solver):" << endl;
	for (int i = 0; i < 3; i++) {
		for (int col = 0; col < 4; col++) {
			if (i == 0) {
				cout << hex << setw(9) << state[i + col * 3] << " ";
			}
			else if (i == 1) {
				cout << hex << setw(9) << state[i + col * 3] << " ";
			}
			else {
				cout << hex << setw(9) << state[i + col * 3] << " ";
			}
		}
		cout << endl;
	}
	cout << endl;
	return true;
}

void Gimli::updateConditionMatrix(UINT32* input, UINT32* output, bool** matrix, int columnNum) {
	bool IX[32], IY[32], IZ[32];
	bool OX[32], OY[32], OZ[32];
	int pos[2], coef[2];
	bool value = 0;
	bool equation[193];
	memset(equation, 0, 193);

	wordToBit(input[0], IX);
	wordToBit(input[1], IY);
	wordToBit(input[2], IZ);
	wordToBit(output[0], OX);
	wordToBit(output[1], OY);
	wordToBit(output[2], OZ);

	int returnedValue = 0;
	for (int depth = 0; depth < 32; depth++) {
		extraLinearRelationsFromX(depth, IX, IY, IZ, OX, pos, coef, value);
		pos[0] = pos[0] + 96 * columnNum;
		pos[1] = pos[1] + 96 * columnNum;
		returnedValue = addRelationToSystem(pos, coef, value, matrix, 192, equation, 193);
		if (returnedValue == 0) {
			cout << "wrong!" << endl;
		}

		extraLinearRelationsFromY(depth, IX, IY, IZ, OY, pos, coef, value);
		pos[0] = pos[0] + 96 * columnNum;
		pos[1] = pos[1] + 96 * columnNum;
		returnedValue = addRelationToSystem(pos, coef, value, matrix, 192, equation, 193);
		if (returnedValue == 0) {
			cout << "wrong!" << endl;
		}

		extraLinearRelationsFromZ(depth, IX, IY, IZ, OZ, pos, coef, value);
		pos[0] = pos[0] + 96 * columnNum;
		pos[1] = pos[1] + 96 * columnNum;
		returnedValue = addRelationToSystem(pos, coef, value, matrix, 192, equation, 193);
		if (returnedValue == 0) {
			cout << "wrong!" << endl;
		}
	}
}

//6-round SFS collision
void Gimli::search6RSFSCollision(int threadNum) {
	GRBEnv env = GRBEnv();
	//env.set(GRB_IntParam_OutputFlag, 0);
	env.set(GRB_IntParam_Threads, threadNum);
	GRBModel model = GRBModel(env);
	vector<vector<GRBVar> > ac;//before constant addition
	vector<vector<vector<GRBVar> > > d;//difference
	vector<vector<vector<GRBVar> > > dn;//difference of nonlinear part
	vector<vector<vector<GRBVar> > > v;//value
	int rounds = 5;

	d.resize(rounds);
	dn.resize(rounds);
	v.resize(rounds);
	ac.resize(6);//6 in total
	int variableNum = 32;

	for (int i = 0; i < 6; i++) {
		ac[i].resize(variableNum);
		for (int j = 0; j < 32; j++) {
			ac[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		}
	}

	for (int i = 0; i < rounds; i++) {
		d[i].resize(12);
		v[i].resize(12);
		dn[i].resize(12);
	}
	for (int r = 0; r < rounds; r++) {
		for (int b = 0; b < 12; b++) {
			d[r][b].resize(variableNum);
			v[r][b].resize(variableNum);
			dn[r][b].resize(variableNum);
		}
	}

	for (int r = 0; r < rounds; r++) {
		for (int b = 0; b < 12; b++) {
			for (int z = 0; z < variableNum; z++) {
				dn[r][b][z] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
				d[r][b][z] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
				v[r][b][z] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
			}
		}
	}

	int r = 0;
	for (int i = 24; i > 20; i--) {
		r = 24 - i;
		if (i % 4 == 0) {//Small swap
			modelXUpdate(model, d[r][0], d[r][1], d[r][2], d[r + 1][3], dn[r][0], v[r][0], v[r][1]);//small swap
			modelYUpdate(model, d[r][0], d[r][1], d[r][2], d[r + 1][1], dn[r][1], v[r][0], v[r][2]);
			modelZUpdate(model, d[r][0], d[r][1], d[r][2], d[r + 1][2], dn[r][2], v[r][1], v[r][2]);

			modelXUpdate(model, d[r][3], d[r][4], d[r][5], d[r + 1][0], dn[r][3], v[r][3], v[r][4]);
			modelYUpdate(model, d[r][3], d[r][4], d[r][5], d[r + 1][4], dn[r][4], v[r][3], v[r][5]);
			modelZUpdate(model, d[r][3], d[r][4], d[r][5], d[r + 1][5], dn[r][5], v[r][4], v[r][5]);

			modelXUpdate(model, d[r][6], d[r][7], d[r][8], d[r + 1][9], dn[r][6], v[r][6], v[r][7]);//small swap
			modelYUpdate(model, d[r][6], d[r][7], d[r][8], d[r + 1][7], dn[r][7], v[r][6], v[r][8]);
			modelZUpdate(model, d[r][6], d[r][7], d[r][8], d[r + 1][8], dn[r][8], v[r][7], v[r][8]);

			modelXUpdate(model, d[r][9], d[r][10], d[r][11], d[r + 1][6], dn[r][9], v[r][9], v[r][10]);
			modelYUpdate(model, d[r][9], d[r][10], d[r][11], d[r + 1][10], dn[r][10], v[r][9], v[r][11]);
			modelZUpdate(model, d[r][9], d[r][10], d[r][11], d[r + 1][11], dn[r][11], v[r][10], v[r][11]);
		}
		else if (i % 2 == 0) {//Big swap
			modelXUpdate(model, d[r][0], d[r][1], d[r][2], d[r + 1][6], dn[r][0], v[r][0], v[r][1]);//big swap
			modelYUpdate(model, d[r][0], d[r][1], d[r][2], d[r + 1][1], dn[r][1], v[r][0], v[r][2]);
			modelZUpdate(model, d[r][0], d[r][1], d[r][2], d[r + 1][2], dn[r][2], v[r][1], v[r][2]);

			modelXUpdate(model, d[r][3], d[r][4], d[r][5], d[r + 1][9], dn[r][3], v[r][3], v[r][4]);
			modelYUpdate(model, d[r][3], d[r][4], d[r][5], d[r + 1][4], dn[r][4], v[r][3], v[r][5]);
			modelZUpdate(model, d[r][3], d[r][4], d[r][5], d[r + 1][5], dn[r][5], v[r][4], v[r][5]);

			modelXUpdate(model, d[r][6], d[r][7], d[r][8], d[r + 1][0], dn[r][6], v[r][6], v[r][7]);//big swap
			modelYUpdate(model, d[r][6], d[r][7], d[r][8], d[r + 1][7], dn[r][7], v[r][6], v[r][8]);
			modelZUpdate(model, d[r][6], d[r][7], d[r][8], d[r + 1][8], dn[r][8], v[r][7], v[r][8]);

			modelXUpdate(model, d[r][9], d[r][10], d[r][11], d[r + 1][3], dn[r][9], v[r][9], v[r][10]);
			modelYUpdate(model, d[r][9], d[r][10], d[r][11], d[r + 1][10], dn[r][10], v[r][9], v[r][11]);
			modelZUpdate(model, d[r][9], d[r][10], d[r][11], d[r + 1][11], dn[r][11], v[r][10], v[r][11]);
		}
		else {//No swap
			modelXUpdate(model, d[r][0], d[r][1], d[r][2], d[r + 1][0], dn[r][0], v[r][0], v[r][1]);//big swap
			modelYUpdate(model, d[r][0], d[r][1], d[r][2], d[r + 1][1], dn[r][1], v[r][0], v[r][2]);
			modelZUpdate(model, d[r][0], d[r][1], d[r][2], d[r + 1][2], dn[r][2], v[r][1], v[r][2]);

			modelXUpdate(model, d[r][3], d[r][4], d[r][5], d[r + 1][3], dn[r][3], v[r][3], v[r][4]);
			modelYUpdate(model, d[r][3], d[r][4], d[r][5], d[r + 1][4], dn[r][4], v[r][3], v[r][5]);
			modelZUpdate(model, d[r][3], d[r][4], d[r][5], d[r + 1][5], dn[r][5], v[r][4], v[r][5]);

			modelXUpdate(model, d[r][6], d[r][7], d[r][8], d[r + 1][6], dn[r][6], v[r][6], v[r][7]);//big swap
			modelYUpdate(model, d[r][6], d[r][7], d[r][8], d[r + 1][7], dn[r][7], v[r][6], v[r][8]);
			modelZUpdate(model, d[r][6], d[r][7], d[r][8], d[r + 1][8], dn[r][8], v[r][7], v[r][8]);

			modelXUpdate(model, d[r][9], d[r][10], d[r][11], d[r + 1][9], dn[r][9], v[r][9], v[r][10]);
			modelYUpdate(model, d[r][9], d[r][10], d[r][11], d[r + 1][10], dn[r][10], v[r][9], v[r][11]);
			modelZUpdate(model, d[r][9], d[r][10], d[r][11], d[r + 1][11], dn[r][11], v[r][10], v[r][11]);
		}

		//model value transition
		if (i % 4 == 0) {//Small swap (and constant addition)
			modelXValueUpdate(model, v[r][0], v[r][1], v[r][2], v[r + 1][3]);//small_swap
			modelYValueUpdate(model, v[r][0], v[r][1], v[r][2], v[r + 1][1]);
			modelZValueUpdate(model, v[r][0], v[r][1], v[r][2], v[r + 1][2]);

			//before v[r+1][0], there is a constant addition and we need to consider it.
			//ac[6-i/4]-> v[r+1][0]
			//modelXValueUpdate(model, v[r][3], v[r][4], v[r][5], v[r + 1][0]);//small_swap
			modelXValueUpdate(model, v[r][3], v[r][4], v[r][5], ac[6 - i / 4]);//small_swap
			//ac -> v[r+1][0]
			for (int j = 0; j < 32; j++) {
				if (constant[6 - i / 4][j] == 0) {
					model.addConstr(ac[6 - i / 4][j] == v[r + 1][0][j]);
				}
				else {
					model.addConstr(ac[6 - i / 4][j] + v[r + 1][0][j] == 1);
				}
			}
			modelYValueUpdate(model, v[r][3], v[r][4], v[r][5], v[r + 1][4]);
			modelZValueUpdate(model, v[r][3], v[r][4], v[r][5], v[r + 1][5]);

			modelXValueUpdate(model, v[r][6], v[r][7], v[r][8], v[r + 1][9]);//small_swap
			modelYValueUpdate(model, v[r][6], v[r][7], v[r][8], v[r + 1][7]);
			modelZValueUpdate(model, v[r][6], v[r][7], v[r][8], v[r + 1][8]);

			modelXValueUpdate(model, v[r][9], v[r][10], v[r][11], v[r + 1][6]);//small_swap
			modelYValueUpdate(model, v[r][9], v[r][10], v[r][11], v[r + 1][10]);
			modelZValueUpdate(model, v[r][9], v[r][10], v[r][11], v[r + 1][11]);
		}
		else if (i % 2 == 0) {
			modelXValueUpdate(model, v[r][0], v[r][1], v[r][2], v[r + 1][6]);//big_swap
			modelYValueUpdate(model, v[r][0], v[r][1], v[r][2], v[r + 1][1]);
			modelZValueUpdate(model, v[r][0], v[r][1], v[r][2], v[r + 1][2]);

			modelXValueUpdate(model, v[r][3], v[r][4], v[r][5], v[r + 1][9]);//big_swap
			modelYValueUpdate(model, v[r][3], v[r][4], v[r][5], v[r + 1][4]);
			modelZValueUpdate(model, v[r][3], v[r][4], v[r][5], v[r + 1][5]);

			modelXValueUpdate(model, v[r][6], v[r][7], v[r][8], v[r + 1][0]);//big_swap
			modelYValueUpdate(model, v[r][6], v[r][7], v[r][8], v[r + 1][7]);
			modelZValueUpdate(model, v[r][6], v[r][7], v[r][8], v[r + 1][8]);

			modelXValueUpdate(model, v[r][9], v[r][10], v[r][11], v[r + 1][3]);//big_swap
			modelYValueUpdate(model, v[r][9], v[r][10], v[r][11], v[r + 1][10]);
			modelZValueUpdate(model, v[r][9], v[r][10], v[r][11], v[r + 1][11]);
		}
		else {
			modelXValueUpdate(model, v[r][0], v[r][1], v[r][2], v[r + 1][0]);
			modelYValueUpdate(model, v[r][0], v[r][1], v[r][2], v[r + 1][1]);
			modelZValueUpdate(model, v[r][0], v[r][1], v[r][2], v[r + 1][2]);

			modelXValueUpdate(model, v[r][3], v[r][4], v[r][5], v[r + 1][3]);
			modelYValueUpdate(model, v[r][3], v[r][4], v[r][5], v[r + 1][4]);
			modelZValueUpdate(model, v[r][3], v[r][4], v[r][5], v[r + 1][5]);

			modelXValueUpdate(model, v[r][6], v[r][7], v[r][8], v[r + 1][6]);
			modelYValueUpdate(model, v[r][6], v[r][7], v[r][8], v[r + 1][7]);
			modelZValueUpdate(model, v[r][6], v[r][7], v[r][8], v[r + 1][8]);

			modelXValueUpdate(model, v[r][9], v[r][10], v[r][11], v[r + 1][9]);
			modelYValueUpdate(model, v[r][9], v[r][10], v[r][11], v[r + 1][10]);
			modelZValueUpdate(model, v[r][9], v[r][10], v[r][11], v[r + 1][11]);
		}
	}

	//input state(S^0)
	for (int i = 1; i < 3; i++) {
		for (int j = 0; j < 4; j++) {
			loadConstraintOnTheOutputSingleDifference(model, d[0][i + j * 3], 0);
		}
	}
	loadConstraintOnTheOutputSingleDifference(model, d[0][0], 0);
	loadConstraintOnTheOutputSingleDifference(model, d[0][6], 0);
	////Test the correctness of the model
	//loadConstraintOnTheOutputSingleDifference(model, d[0][3], 0x7c2c642a);
	////

	//S^1
	loadConstraintOnTheOutputSingleDifference(model, d[1][0], 0);
	loadConstraintOnTheOutputSingleDifference(model, d[1][6], 0);

	//S^0, S^1, S^2, S^3, S^4
	for (int k = 0; k < 5; k++) {
		for (int j = 0; j < 3; j++) {
			for (int z = 0; z < 32; z++) {
				model.addConstr(d[k][j + 3][z] == d[k][j + 9][z]);
			}
		}
	}

	//additional conditions on S^4
	loadConstraintOnTheOutputSingleDifference(model, d[4][3], 0x80);
	loadConstraintOnTheOutputSingleDifference(model, d[4][4], 0x400000);
	loadConstraintOnTheOutputSingleDifference(model, d[4][5], 0x80000000);
	loadConstraintOnTheOutputSingleDifference(model, d[4][9], 0x80);
	loadConstraintOnTheOutputSingleDifference(model, d[4][10], 0x400000);
	loadConstraintOnTheOutputSingleDifference(model, d[4][11], 0x80000000);

	//additional constraints on S^3
	GRBLinExpr sum3 = 0;
	for (int j = 0; j < 32; j++) {
		sum3 += d[3][3][j];
		sum3 += d[3][4][j];
		sum3 += d[3][5][j];
	}
	model.addConstr(sum3 <= 8);

	//correctness test (solution)
	/*loadConstraintOnTheOutputSingleDifference(model, d[0][3], 0x7c2c642a);

	loadConstraintOnTheOutputSingleDifference(model, d[1][4], 0x6e1c342c);
	loadConstraintOnTheOutputSingleDifference(model, d[1][5], 0x2a7c2c64);

	loadConstraintOnTheOutputSingleDifference(model, d[2][3], 0x91143078);
	loadConstraintOnTheOutputSingleDifference(model, d[2][4], 0x28785014);
	loadConstraintOnTheOutputSingleDifference(model, d[2][5], 0x35288a58);

	loadConstraintOnTheOutputSingleDifference(model, d[3][3], 0x80010008);
	loadConstraintOnTheOutputSingleDifference(model, d[3][4], 0x00002000);
	loadConstraintOnTheOutputSingleDifference(model, d[3][5], 0x44400080);*/

	model.optimize();

	UINT32 value[12];
	UINT32 state[12];
	cout << "Discovered trail:" << endl;
	for (int r = 0; r < rounds; r++) {
		for (int b = 0; b < 12; b++) {
			value[b] = 0;
			for (int i = 0; i < 32; i++) {
				if (d[r][b][i].get(GRB_DoubleAttr_X) == 1) {
					value[b] = value[b] | EXP[i];
				}
			}
		}
		//output value
		for (int i = 0; i < 3; i++) {
			for (int col = 0; col < 4; col++) {
				cout << hex << setw(9) << value[i + col * 3] << " ";
			}
			cout << endl;
		}
		cout << endl;
	}
	cout << endl;

	cout << "The message satisfying the differential conditions:" << endl;
	for (int r = 0; r < rounds; r++) {
		for (int b = 0; b < 12; b++) {
			value[b] = 0;
			for (int i = 0; i < 32; i++) {
				if (v[r][b][i].get(GRB_DoubleAttr_X) == 1) {
					value[b] = value[b] | EXP[i];
				}
			}
		}
		if (r == 0) {
			for (int i = 0; i < 12; i++) {
				state[i] = value[i];
			}
		}
		//output value
		for (int i = 0; i < 12; i++) {
			cout << hex << value[i] << " ";
		}
		cout << endl;
	}
	cout << endl;
}

UINT32 Gimli::toUINT32(bool* arr) {
	UINT32 res = 0;
	for (int i = 0; i < 32; i++) {
		if (arr[i]) {
			res |= EXP[i];
		}
	}
	return res;
}

void Gimli::generateValidCapacityFor6RCollision() {
	UINT32 inputDiff[5][3] = {
		{0x7c2c642a,0,0},
		{0,0x6e1c342c,0x2a7c2c64},
		{0x91143078,0x28785014,0x35288a58},
		{0x80010008,0x2000,0x44400080},
		{0x80,0x400000,0x80000000}
	};

	bool** matrix0;
	matrix0 = new bool* [192];
	for (int i = 0; i < 192; i++) {
		matrix0[i] = new bool[193];
	}
	clearMatrix(matrix0, 192, 193);

	UINT32 input[3], output[3];

	input[0] = inputDiff[0][0];
	input[1] = inputDiff[0][1];
	input[2] = inputDiff[0][2];
	output[0] = inputDiff[1][0];
	output[1] = inputDiff[1][1];
	output[2] = inputDiff[1][2];
	updateConditionMatrix(input, output, matrix0, 0);//The condition on S^0
	outputMatrix(matrix0, 192, 193);
	vector<int> freeBitPos;
	freeBitPos.clear();
	bool bitValue[64];
	for (int i = 64; i < 96; i++) {
		if (matrix0[i][i] == 0) {
			freeBitPos.push_back(i-64);
		}
		else {
			bitValue[i - 64] = matrix0[i][192];
		}
	}

	ifstream inFile("S1_Sol.txt");
	Column *validSol;
	int size = 0xcc * 8;
	validSol = new Column[size];
	int counter = 0;
	int diff = 0;
	while (inFile >> validSol[counter].x) {
		if (counter >= size) {
			cout << "wrong" << endl;
			break;
		}
		inFile >> validSol[counter].y >> validSol[counter].z;
		counter++;
	}
	sort(validSol, validSol + counter, compare);

	int count = 0;
	bool flag = false;
	vector<Column> filteredSol;
	filteredSol.clear();
	for (int i = 0; i < counter-1; i++) {
		count++;
		if (validSol[i].y == validSol[i + 1].y && validSol[i].z == validSol[i + 1].z) {
			flag = true;
		}
		else {
			cout << hex << validSol[i].x <<" "<< validSol[i].y << " " << validSol[i].z << ": ";
			cout << count << endl;
			count = 0;
			diff++;
			filteredSol.push_back(validSol[i]);
			flag = false;
		}
	}
	if (flag) {
		cout << hex << validSol[counter-1].x << " " << validSol[counter - 1].y << " " << validSol[counter - 1].z << ": ";
		cout << count << endl;
		filteredSol.push_back(validSol[counter - 1]);
		diff++;
	}
	else {
		cout << hex << validSol[counter - 1].x << " " << validSol[counter - 1].y << " " << validSol[counter - 1].z << ": ";
		cout << 1 << endl;
		filteredSol.push_back(validSol[counter - 1]);
		diff++;
	}
	cout << dec << "diff: " << diff << endl;//0x2d0
	cout << dec << "diff: " << filteredSol.size() << endl;

	//exhaust the capacity part
	int frSize = freeBitPos.size();
	UINT64 end = EXP[frSize];
	UINT32 validNum = 0;

	vector<vector<bool> > YZPrimeSum,ZPrimeBit;
	YZPrimeSum.resize(diff);
	ZPrimeBit.resize(diff);
	for (int i = 0; i < diff; i++) {
		for (int j = 0; j < 32; j++) {
			YZPrimeSum[i].push_back(((filteredSol[i].y ^ filteredSol[i].z) >> j) & 0x1);
			ZPrimeBit[i].push_back((filteredSol[i].z >> j) & 0x1);
		}
	}

	bool equation[193];
	memset(equation, 0, 193);
	flag = true;
	bool **tmpMatrix0;
	tmpMatrix0 = new bool* [192];
	for (int i = 0; i < 192; i++) {
		tmpMatrix0[i] = new bool[193];
	}
	matrixEqual(matrix0, tmpMatrix0, 192, 193);

	cout << hex<< "exhausting space size: 0x" << end << endl;

	vector<Column> S0ValidCap;
	Column validValue;
	for (UINT64 solution = 0; solution < end; solution++) {
		for (int j = 0; j < frSize; j++) {
			bitValue[freeBitPos[j]] = (solution >> j) & 0x1;
		}
		//matrixEqual(tmpMatrix0, matrix0, 192, 193);

		//cout << hex << solution << endl;
		if (solution % EXP[10] == 0) {
			cout << "currectly exhausted size: 0x" << solution<< endl;
		}

		for (int k = 0; k < diff; k++) {
			//matrixEqual(tmpMatrix0, matrix0, 192, 193);
			for (int row = 32; row < 64; row++) {
				for (int col = 32; col < 64; col++) {
					matrix0[row][col] = tmpMatrix0[row][col];
				}
				matrix0[row][192] = tmpMatrix0[row][192];
			}
			
			flag = true;

			//YZPrimeSum[0]
			equation[192] = YZPrimeSum[k][0];
			equation[YMap(MOD(0-9))] = 1;
			flag=addToEquationSystem(equation, matrix0, 192, 193);
			memset(equation, 0, 193);
			if (!flag) {
				continue;
			}

			//YZPrimeSum[1]
			equation[192] = YZPrimeSum[k][1]^(ZPrimeBit[k][0]&(~bitValue[0]));
			equation[YMap(MOD(1 - 9))] = 1;
			flag = addToEquationSystem(equation, matrix0, 192, 193);
			memset(equation, 0, 193);
			if (!flag) {
				continue;
			}

			//YZPrimeSum[2]
			equation[192] = YZPrimeSum[k][2] ^ ((~bitValue[1]) & (ZPrimeBit[k][1] ^ bitValue[0]));
			equation[YMap(MOD(2 - 9))] = 1;
			if(bitValue[0])
				equation[YMap(MOD(0 - 9))] = 1;
			flag = addToEquationSystem(equation, matrix0, 192, 193);
			memset(equation, 0, 193);
			if (!flag) {
				continue;
			}

			//YZPrimeSum[3-31]
			for (int d = 3; d < 32; d++) {
				equation[192] = YZPrimeSum[k][d];
				equation[YMap(MOD(d - 9))] = 1;
				if (bitValue[d - 2]) {
					equation[YMap(MOD(d - 11))] = 1;
				}
				if (bitValue[d - 1] == 0) {
					equation[192] = equation[192] ^ ZPrimeBit[k][d - 1] ^ bitValue[d - 2];
					if (bitValue[d - 3]) {
						equation[YMap(MOD(d - 12))] = 1;
					}
				}
				flag=addToEquationSystem(equation, matrix0, 192, 193);
				memset(equation, 0, 193);
				if (flag == false) {
					break;
				}
			}

			if (flag) {
				//outputMatrix(matrix0, 192, 193);
				//system("pause");
				//counter the number of solutions for y
				UINT32 solNum = 0;
				UINT32 IY=0, IZ=0;
				IZ = toUINT32(bitValue);
				for (int m = 32; m < 64; m++) {
					if (matrix0[m][m] == 0) {
						solNum++;
					}
					else {
						if (matrix0[m][192]) {
							IY = IY | EXP[m - 32];
						}
					}
				}
				validValue.x = 0;
				validValue.y = IY;
				validValue.z = IZ;
				S0ValidCap.push_back(validValue);
				//IY,IZ,OY,OZ are known, we can check the correctness
				/*IY = LL(IY, 9);
				UINT32 IX, OX, OY, OZ;
				OY = filteredSol[k].y;
				OZ = filteredSol[k].z;
				IX = OZ ^ (IZ << 1) ^ ((IY & IZ) << 2);
				OY = IY ^ IX ^ ((IX | IZ) << 1);

				validNum += EXP[solNum];
				//validNum++;
				cout << "correct: " << solNum << " " << IY << " " << IZ;
				if (OY == filteredSol[k].y) {
					cout << " right" << endl;
				}
				else {
					cout << " wrong" << endl;
				}*/
			}
		}
	}
	cout << "validNum: 0x" << validNum << endl;
	cout << "validNum: 0x" << S0ValidCap.size() << endl;
	//write to file
	ofstream outS0("S0.txt");
	for (int i = 0; i < S0ValidCap.size(); i++) {
		outS0 << dec << S0ValidCap[i].y << " " << S0ValidCap[i].z << endl;
	}
	outS0.close();

	for (int i = 0; i < 192; i++) {
		delete[]matrix0[i];
		delete[]tmpMatrix0[i];
	}
	delete[]matrix0;
	delete[]tmpMatrix0;
	freeBitPos.clear();
}

void Gimli::countValidCapacityFor6RCollision() {
	ifstream inFile("S0.txt");
	Column* validSol;
	int size = 0x6990;
	validSol = new Column[size];
	int diff = 0;
	for (int i = 0; i < size; i++) {
		inFile >> validSol[i].y >> validSol[i].z;
		validSol[i].x = 0;
	}
	inFile.close();
	sort(validSol, validSol + size, compare);

	int count = 0;
	bool flag = false;
	vector<Column> filteredSol;
	filteredSol.clear();
	for (int i = 0; i < size - 1; i++) {
		count++;
		if (validSol[i].y == validSol[i + 1].y && validSol[i].z == validSol[i + 1].z) {
			flag = true;
		}
		else {
			cout << hex << validSol[i].y << " " << validSol[i].z << ": ";
			cout << count << endl;
			count = 0;
			diff++;
			filteredSol.push_back(validSol[i]);
			flag = false;
		}
	}
	if (flag) {
		cout << hex << " " << validSol[size - 1].y << " " << validSol[size - 1].z << ": ";
		cout << 2 << endl;
		filteredSol.push_back(validSol[size - 1]);
		diff++;
	}
	else {
		cout << " " << validSol[size - 2].y << " " << validSol[size - 2].z << endl;
		cout << " " << validSol[size - 1].y << " " << validSol[size - 1].z << endl;
		filteredSol.push_back(validSol[size - 1]);
		diff++;
	}
	cout << hex << "diff: " << diff << endl;//0x34c8
	cout << hex << "diff: " << filteredSol.size() << endl;

	delete[]validSol;
	filteredSol.clear();
}

//14-round zero-internal-differential distinguisher
bool Gimli::zeroInternalDiffAttack(UINT32 constant[]) {
	//First consider the following sequence of operations (9 rounds)
	//(SP)->(SP->B_SW)->(SP)->(SP->S_SW->AC)->
	//(SP)->(SP->B_SW)->(SP)->(SP->S_SW->AC)->
	//(SP)
	bool zeroDiffBit[6];

	UINT32 state[12],tmpState[12];
	for (int i = 0; i < 3; i++) {
		state[i] = getRand32();
		state[i + 3] = state[i];
		state[i + 6] = state[i];
		state[i + 9] = state[i];
	}
	for (int i = 0; i < 12; i++) {
		tmpState[i] = state[i];
	}
	//Apply the permutation
	for (int i = 0; i < 2; i++) {
		SPColumn(state, 12);

		SPColumn(state, 12);
		bigSwap(state);

		SPColumn(state, 12);

		SPColumn(state, 12);
		smallSwap(state);
		state[0] ^= constant[i];
	}
	SPColumn(state, 12);
	//compute the 4 zero-internal-difference bits
	//(state[3]^state[9])[0,1,2]=0
	//state[4][0]^state[5][0]=state[10][0]^state[11][0]
	zeroDiffBit[0] = BIT(state[3], 0) ^ BIT(state[9], 0);
	zeroDiffBit[1] = BIT(state[3], 1) ^ BIT(state[9], 1);
	zeroDiffBit[2] = BIT(state[3], 2) ^ BIT(state[9], 2);
	zeroDiffBit[3] = BIT(state[4], 0) ^ BIT(state[5], 0) ^ BIT(state[10], 0) ^ BIT(state[11], 0);

	//Next consider the following sequence of operations
	//(SP->S_SW->AC)->(SP)->(SP->B_SW)->(SP)->(SP->S_SW->AC) (5 rounds)
	tmpState[0] ^= constant[2];
	smallSwap(tmpState);
	inverseSPColumn(tmpState, 12);
	
	inverseSPColumn(tmpState, 12);

	bigSwap(tmpState);
	inverseSPColumn(tmpState, 12);

	inverseSPColumn(tmpState, 12);

	tmpState[0] ^= constant[3];
	smallSwap(tmpState);
	inverseSPColumn(tmpState, 12);

	//distinguisher: tmpState[0][8]=tmpState[6][8]
	//distinguisher: tmpState[1][23]=tmpState[7][23]
	zeroDiffBit[4] = BIT(tmpState[0], 8) ^ BIT(tmpState[6], 8);
	zeroDiffBit[5] = BIT(tmpState[1], 23) ^ BIT(tmpState[7], 23);

	/*for (int i = 1; i < 7; i++) {
		cout << "distinguisher " << i<<": " << zeroDiffBit[i-1] << endl;
	}*/

	for (int i = 0; i < 6; i++) {
		if (zeroDiffBit[i]) {
			return false;
		}
	}
	return true;
}

void Gimli::verifySolForS1() {
	ifstream inFile("S1_Sol.txt");
	Column* validSol;
	int size = 0xcc * 8;
	validSol = new Column[size];
	int counter = 0;
	while (inFile >> validSol[counter].x) {
		if (counter >= size) {
			cout << "wrong" << endl;
			break;
		}
		inFile >> validSol[counter].y >> validSol[counter].z;
		counter++;
	}
	inFile.close();

	UINT32 inputDiff[5][3] = {
		{0x7c2c642a,0,0},
		{0,0x6e1c342c,0x2a7c2c64},
		{0x91143078,0x28785014,0x35288a58},
		{0x80010008,0x2000,0x44400080},
		{0x80,0x400000,0x80000000}
	};
	
	//UINT32 state[3],stateDiff[3];
	UINT32 state[6], stateDiff[6],tmp;
	bool flag = true;
	for (int i = 0; i < counter; i++) {
		for (int j = 0; j < counter; j++) {
			state[0] = validSol[i].x;
			state[1] = validSol[i].y;
			state[2] = validSol[i].z;
			stateDiff[0] = state[0] ^ inputDiff[1][0];
			stateDiff[1] = state[1] ^ inputDiff[1][1];
			stateDiff[2] = state[2] ^ inputDiff[1][2];

			state[3]= validSol[j].x;
			state[4] = validSol[j].y;
			state[5] = validSol[j].z;
			stateDiff[3] = state[3] ^ inputDiff[1][0];
			stateDiff[4] = state[4] ^ inputDiff[1][1];
			stateDiff[5] = state[5] ^ inputDiff[1][2];

			SPColumn(state, 6);//SP
			SPColumn(state, 6);//SP->B_Swap
			tmp = state[0];
			state[0] = state[3];
			state[3] = state[0];
			SPColumn(state, 6);//SP

			SPColumn(stateDiff, 6);//SP
			SPColumn(stateDiff, 6);//SP->B_Swap
			tmp = stateDiff[0];
			stateDiff[0] = stateDiff[3];
			stateDiff[3] = stateDiff[0];
			SPColumn(stateDiff, 6);//SP

			for (int j = 0; j < 6; j++) {
				if ((state[j] ^ stateDiff[j]) != inputDiff[4][j % 3]) {
					cout << "wrong" << endl;
					flag = false;
				}
			}
		}
	}
	if (flag) {
		cout << "All are correct!" << endl;
	}
}

//Check the validity of the differential pattern (8 rounds)
void Gimli::searchCollision4RPattern() {
	GRBEnv env = GRBEnv();
	//env.set(GRB_IntParam_OutputFlag, 0);
	env.set(GRB_IntParam_Threads, 3);
	GRBModel model = GRBModel(env);
	vector<vector<vector<GRBVar> > > d;//difference
	vector<vector<vector<GRBVar> > > dn;//difference of nonlinear part
	vector<vector<vector<GRBVar> > > v;//value
	int rounds = 5;

	d.resize(rounds);
	dn.resize(rounds);
	v.resize(rounds);

	int variableNum = 32;
	for (int i = 0; i < rounds; i++) {
		d[i].resize(6);
		v[i].resize(6);
		dn[i].resize(6);
	}

	for (int r = 0; r < rounds; r++) {
		for (int b = 0; b < 6; b++) {
			d[r][b].resize(variableNum);
			v[r][b].resize(variableNum);
			dn[r][b].resize(variableNum);
		}
	}

	for (int r = 0; r < rounds; r++) {
		for (int b = 0; b < 6; b++) {
			for (int z = 0; z < variableNum; z++) {
				dn[r][b][z] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
				d[r][b][z] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
				v[r][b][z] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
			}
		}
	}

	//write inequalities
	//SP->(SP-B_SW)->(SP)->(SP)
	for (int r = 0; r < 4; r++) {
		if (r == 0 || r==2 || r==3) {
			modelXUpdate(model, d[r][0], d[r][1], d[r][2], d[r + 1][0], dn[r][0], v[r][0], v[r][1]);
			modelYUpdate(model, d[r][0], d[r][1], d[r][2], d[r + 1][1], dn[r][1], v[r][0], v[r][2]);
			modelZUpdate(model, d[r][0], d[r][1], d[r][2], d[r + 1][2], dn[r][2], v[r][1], v[r][2]);

			modelXValueUpdate(model, v[r][0], v[r][1], v[r][2], v[r + 1][0]);
			modelYValueUpdate(model, v[r][0], v[r][1], v[r][2], v[r + 1][1]);
			modelZValueUpdate(model, v[r][0], v[r][1], v[r][2], v[r + 1][2]);

			modelXUpdate(model, d[r][3], d[r][4], d[r][5], d[r + 1][3], dn[r][3], v[r][3], v[r][4]);
			modelYUpdate(model, d[r][3], d[r][4], d[r][5], d[r + 1][4], dn[r][4], v[r][3], v[r][5]);
			modelZUpdate(model, d[r][3], d[r][4], d[r][5], d[r + 1][5], dn[r][5], v[r][4], v[r][5]);

			modelXValueUpdate(model, v[r][3], v[r][4], v[r][5], v[r + 1][3]);
			modelYValueUpdate(model, v[r][3], v[r][4], v[r][5], v[r + 1][4]);
			modelZValueUpdate(model, v[r][3], v[r][4], v[r][5], v[r + 1][5]);
		}
		if (r == 1) {//Big-Swap
			modelXUpdate(model, d[r][0], d[r][1], d[r][2], d[r + 1][3], dn[r][0], v[r][0], v[r][1]);
			modelYUpdate(model, d[r][0], d[r][1], d[r][2], d[r + 1][1], dn[r][1], v[r][0], v[r][2]);
			modelZUpdate(model, d[r][0], d[r][1], d[r][2], d[r + 1][2], dn[r][2], v[r][1], v[r][2]);

			modelXValueUpdate(model, v[r][0], v[r][1], v[r][2], v[r + 1][3]);
			modelYValueUpdate(model, v[r][0], v[r][1], v[r][2], v[r + 1][1]);
			modelZValueUpdate(model, v[r][0], v[r][1], v[r][2], v[r + 1][2]);

			modelXUpdate(model, d[r][3], d[r][4], d[r][5], d[r + 1][0], dn[r][3], v[r][3], v[r][4]);
			modelYUpdate(model, d[r][3], d[r][4], d[r][5], d[r + 1][4], dn[r][4], v[r][3], v[r][5]);
			modelZUpdate(model, d[r][3], d[r][4], d[r][5], d[r + 1][5], dn[r][5], v[r][4], v[r][5]);

			modelXValueUpdate(model, v[r][3], v[r][4], v[r][5], v[r + 1][0]);
			modelYValueUpdate(model, v[r][3], v[r][4], v[r][5], v[r + 1][4]);
			modelZValueUpdate(model, v[r][3], v[r][4], v[r][5], v[r + 1][5]);
		}
	}

	GRBLinExpr sum0 = 0;
	for (int i = 0; i < 32; i++) {
		sum0 += d[0][3][i];
	}
	model.addConstr(sum0 >= 1);

	loadConstraintOnTheOutputSingleDifference(model, d[0][0], 0);
	loadConstraintOnTheOutputSingleDifference(model, d[0][1], 0);
	loadConstraintOnTheOutputSingleDifference(model, d[0][2], 0);
	loadConstraintOnTheOutputSingleDifference(model, d[0][4], 0);
	loadConstraintOnTheOutputSingleDifference(model, d[0][5], 0);

	loadConstraintOnTheOutputSingleDifference(model, d[2][0], 0);

	loadConstraintOnTheOutputSingleDifference(model, d[4][4], 0);
	loadConstraintOnTheOutputSingleDifference(model, d[4][5], 0);

	model.optimize();

	UINT32 value[6];
	UINT32 state[6];
	for (int r = 0; r < rounds; r++) {
		for (int b = 0; b < 6; b++) {
			value[b] = 0;
			for (int i = 0; i < 32; i++) {
				if (d[r][b][i].get(GRB_DoubleAttr_X) == 1) {
					value[b] = value[b] | EXP[i];
				}
			}
		}
		//output value
		for (int i = 0; i < 6; i++) {
			cout << hex << value[i] << " ";
		}
		cout << endl;
	}
	cout << endl;

	for (int r = 0; r < rounds; r++) {
		for (int b = 0; b < 6; b++) {
			value[b] = 0;
			for (int i = 0; i < 32; i++) {
				if (v[r][b][i].get(GRB_DoubleAttr_X) == 1) {
					value[b] = value[b] | EXP[i];
				}
			}
		}
		//output value
		for (int i = 0; i < 6; i++) {
			cout << hex << value[i] << " ";
		}
		cout << endl;
	}
	cout << endl;
}

//Check the validity of the differential pattern (7 rounds)
void Gimli::searchCollision3RPattern() {
	GRBEnv env = GRBEnv();
	//env.set(GRB_IntParam_OutputFlag, 0);
	env.set(GRB_IntParam_Threads, 3);
	GRBModel model = GRBModel(env);
	vector<vector<vector<GRBVar> > > d;//difference
	vector<vector<vector<GRBVar> > > dn;//difference of nonlinear part
	vector<vector<vector<GRBVar> > > v;//value
	int rounds = 4;

	d.resize(rounds);
	dn.resize(rounds);
	v.resize(rounds);

	int variableNum = 32;
	for (int i = 0; i < rounds; i++) {
		d[i].resize(6);
		v[i].resize(6);
		dn[i].resize(6);
	}

	for (int r = 0; r < rounds; r++) {
		for (int b = 0; b < 6; b++) {
			d[r][b].resize(variableNum);
			v[r][b].resize(variableNum);
			dn[r][b].resize(variableNum);
		}
	}

	for (int r = 0; r < rounds; r++) {
		for (int b = 0; b < 6; b++) {
			for (int z = 0; z < variableNum; z++) {
				dn[r][b][z] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
				d[r][b][z] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
				v[r][b][z] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
			}
		}
	}

	//write inequalities
	//SP->(SP-B_SW)->(SP)
	for (int r = 0; r < 3; r++) {
		if (r == 0 || r == 2) {
			modelXUpdate(model, d[r][0], d[r][1], d[r][2], d[r + 1][0], dn[r][0], v[r][0], v[r][1]);
			modelYUpdate(model, d[r][0], d[r][1], d[r][2], d[r + 1][1], dn[r][1], v[r][0], v[r][2]);
			modelZUpdate(model, d[r][0], d[r][1], d[r][2], d[r + 1][2], dn[r][2], v[r][1], v[r][2]);

			modelXValueUpdate(model, v[r][0], v[r][1], v[r][2], v[r + 1][0]);
			modelYValueUpdate(model, v[r][0], v[r][1], v[r][2], v[r + 1][1]);
			modelZValueUpdate(model, v[r][0], v[r][1], v[r][2], v[r + 1][2]);

			modelXUpdate(model, d[r][3], d[r][4], d[r][5], d[r + 1][3], dn[r][3], v[r][3], v[r][4]);
			modelYUpdate(model, d[r][3], d[r][4], d[r][5], d[r + 1][4], dn[r][4], v[r][3], v[r][5]);
			modelZUpdate(model, d[r][3], d[r][4], d[r][5], d[r + 1][5], dn[r][5], v[r][4], v[r][5]);

			modelXValueUpdate(model, v[r][3], v[r][4], v[r][5], v[r + 1][3]);
			modelYValueUpdate(model, v[r][3], v[r][4], v[r][5], v[r + 1][4]);
			modelZValueUpdate(model, v[r][3], v[r][4], v[r][5], v[r + 1][5]);
		}
		if (r == 1) {//Big-Swap
			modelXUpdate(model, d[r][0], d[r][1], d[r][2], d[r + 1][3], dn[r][0], v[r][0], v[r][1]);
			modelYUpdate(model, d[r][0], d[r][1], d[r][2], d[r + 1][1], dn[r][1], v[r][0], v[r][2]);
			modelZUpdate(model, d[r][0], d[r][1], d[r][2], d[r + 1][2], dn[r][2], v[r][1], v[r][2]);

			modelXValueUpdate(model, v[r][0], v[r][1], v[r][2], v[r + 1][3]);
			modelYValueUpdate(model, v[r][0], v[r][1], v[r][2], v[r + 1][1]);
			modelZValueUpdate(model, v[r][0], v[r][1], v[r][2], v[r + 1][2]);

			modelXUpdate(model, d[r][3], d[r][4], d[r][5], d[r + 1][0], dn[r][3], v[r][3], v[r][4]);
			modelYUpdate(model, d[r][3], d[r][4], d[r][5], d[r + 1][4], dn[r][4], v[r][3], v[r][5]);
			modelZUpdate(model, d[r][3], d[r][4], d[r][5], d[r + 1][5], dn[r][5], v[r][4], v[r][5]);

			modelXValueUpdate(model, v[r][3], v[r][4], v[r][5], v[r + 1][0]);
			modelYValueUpdate(model, v[r][3], v[r][4], v[r][5], v[r + 1][4]);
			modelZValueUpdate(model, v[r][3], v[r][4], v[r][5], v[r + 1][5]);
		}
	}

	GRBLinExpr sum0 = 0;
	for (int i = 0; i < 32; i++) {
		sum0 += d[0][3][i];
	}
	model.addConstr(sum0 >= 1);

	loadConstraintOnTheOutputSingleDifference(model, d[0][0], 0);
	loadConstraintOnTheOutputSingleDifference(model, d[0][1], 0);
	loadConstraintOnTheOutputSingleDifference(model, d[0][2], 0);
	loadConstraintOnTheOutputSingleDifference(model, d[0][4], 0);
	loadConstraintOnTheOutputSingleDifference(model, d[0][5], 0);

	loadConstraintOnTheOutputSingleDifference(model, d[2][0], 0);

	loadConstraintOnTheOutputSingleDifference(model, d[3][4], 0);
	loadConstraintOnTheOutputSingleDifference(model, d[3][5], 0);

	model.optimize();
}
