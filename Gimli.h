#ifndef GIMLI_H_
#define GIMLI_H_
#include <vector>
#include <string>
#include "gurobi_c++.h"
using namespace std;
typedef unsigned int UINT32;
typedef unsigned long long UINT64;
#define SIZE 384

#define BIT(x,i) ((x>>i)&0x1)
#define MOD(i) ((i+32)%32)
#define YMap(p) (p+32)
#define ZMap(p) (p+64)

const UINT32 EXP[32] = {
	0x1,0x2,0x4,0x8,
	0x10,0x20,0x40,0x80,
	0x100,0x200,0x400,0x800,
	0x1000,0x2000,0x4000,0x8000,
	0x10000,0x20000,0x40000,0x80000,
	0x100000,0x200000,0x400000,0x800000,
	0x1000000,0x2000000,0x4000000,0x8000000,
	0x10000000,0x20000000,0x40000000,0x80000000
};

struct Block {
	UINT32 block[6];
};

struct Expression
{
	bool constant;
	int var;
	bool isCon;
};

struct Column
{
	UINT32 x, y, z;
};


class Gimli{
public:
	Gimli();
	~Gimli();
	UINT32 LL(UINT32 number, int n);
	UINT32 RR(UINT32 number, int n);
	UINT64 getRand64();
	UINT32 getRand32();
	void SPBox(UINT32 *input, UINT32 *output);
	void SPColumn(UINT32 state[], int size);
	void inverseSPBox(UINT32 input[]);
	void inverseSPColumn(UINT32 state[], int size);
	void bigSwap(UINT32 state[]);
	void smallSwap(UINT32 state[]);
	void initializeMinimalModel(string filename,int **matrix);

	//modeling
	void modelXUpdate(GRBModel &model, vector<GRBVar>& id, vector<GRBVar>& od, vector<GRBVar>& odn, vector<GRBVar>& v);
	void modelYUpdate(GRBModel &model, vector<GRBVar>& id, vector<GRBVar>& od, vector<GRBVar>& odn, vector<GRBVar>& v);
	void modelZUpdate(GRBModel &model, vector<GRBVar>& id, vector<GRBVar>& od, vector<GRBVar>& odn, vector<GRBVar>& v);

	//reload the above three functions
	void modelXUpdate(GRBModel& model, vector<GRBVar>& idX, vector<GRBVar>& idY, vector<GRBVar>& idZ, vector<GRBVar>& odX, vector<GRBVar>& odnX, vector<GRBVar>& vX, vector<GRBVar>& vY);
	void modelYUpdate(GRBModel& model, vector<GRBVar>& idX, vector<GRBVar>& idY, vector<GRBVar>& idZ, vector<GRBVar>& odY, vector<GRBVar>& odnY, vector<GRBVar>& vX, vector<GRBVar>& vZ);
	void modelZUpdate(GRBModel& model, vector<GRBVar>& idX, vector<GRBVar>& idY, vector<GRBVar>& idZ, vector<GRBVar>& odZ, vector<GRBVar>& odnZ, vector<GRBVar>& vY, vector<GRBVar>& vZ);

	void modelXValueUpdate(GRBModel& model, vector<GRBVar>& iv, vector<GRBVar>& ov);
	void modelYValueUpdate(GRBModel& model, vector<GRBVar>& iv, vector<GRBVar>& ov);
	void modelZValueUpdate(GRBModel& model, vector<GRBVar>& iv, vector<GRBVar>& ov);

	//reload the above three functions
	void modelXValueUpdate(GRBModel& model, vector<GRBVar>& iX, vector<GRBVar>& iY, vector<GRBVar>& iZ, vector<GRBVar>& oX);
	void modelYValueUpdate(GRBModel& model, vector<GRBVar>& iX, vector<GRBVar>& iY, vector<GRBVar>& iZ, vector<GRBVar>& oY);
	void modelZValueUpdate(GRBModel& model, vector<GRBVar>& iX, vector<GRBVar>& iY, vector<GRBVar>& iZ, vector<GRBVar>& oZ);

	void loadConstraintThreeLinearZero(GRBVar& a, GRBVar& b, GRBVar& c, GRBModel& model);
	void loadConstraintFourLinearZero(GRBVar& a, GRBVar& b, GRBVar& c, GRBVar& d, GRBModel& model);
	void loadAndConstraint(GRBVar& a, GRBVar& b, GRBVar& c, GRBVar& d, GRBVar& e, GRBModel& model);
	void loadORConstraint(GRBVar& a, GRBVar& b, GRBVar& c, GRBVar& d, GRBVar& e, GRBModel& model);

	void loadConstraintOnTheOutputDifference(GRBModel& model, vector<GRBVar>& di, UINT32 output[]);
	void loadConstraintOnTheOutputSingleDifference(GRBModel& model, vector<GRBVar>& di, UINT32 output);

	//constraints for value
	void loadAndValueConstraint(GRBVar& a, GRBVar& b, GRBVar& c, GRBVar& d, GRBVar& e, GRBModel& model);
	void loadORValueConstraint(GRBVar& a, GRBVar& b, GRBVar& c, GRBVar& d, GRBVar& e, GRBModel& model);

	//model the whole state
	void findCorrectDifferential(int bound);

	//find compatible message pair using message modification
	void wordToBit(UINT32 word, bool* bitArray);
	void extraLinearRelationsFromX(int bitPos, bool* IX, bool* IY, bool* IZ, bool* OX,int pos[],int coef[],bool &value);
	void extraLinearRelationsFromY(int bitPos, bool* IX, bool* IY, bool* IZ, bool* OY, int pos[], int coef[], bool &value);
	void extraLinearRelationsFromZ(int bitPos, bool* IX, bool* IY, bool* IZ, bool* OZ, int pos[], int coef[], bool &value);
	void outputCondition(int pos[],int coef[],bool value,int columnNum);

	//process constructing the equation system
	int addToEquationSystem(bool equation[], bool** equationSystem, int row, int column);
	void outputMatrix(bool** matrix, int row, int column);
	int countFreedom(bool** matrix, int row, int column);
	void clearMatrix(bool** matrix, int row, int column);
	void matrixEqual(bool** sourceMatrix, bool** desMatrix, int row, int column);
	
	int addRelationToSystem(int pos[], int coef[], bool& value,bool **system, int row, bool *equation,int eqSize);

	void checkRealSFSCollision();
	void permutation(UINT32 state[],int rounds);
	void permutation(UINT32 state[], int start,int end);
	void outputState(UINT32 state[]);

	//Test the correctness of a given diff (e.g. the official 12-round differential)
	bool testCorrectOfDiff(UINT32 inputDiff[][12],int rounds,int start,int end);
	void updateConditionMatrix(UINT32* input, UINT32* output, bool** matrix, int columnNum);

	//search 6-round diff
	void search6RSFSCollision(int threadNum);

	string toBinary(UINT32 num);
	UINT32 toUINT32(bool* arr);
	void generateValidCapacityFor6RCollision();
	void countValidCapacityFor6RCollision();
	void verifySolForS1();//verify validity of the solution for the second column of S^1
	void searchCollision3RPattern();//Check the validity of the differential pattern (7 rounds)
	void searchCollision4RPattern();//Check the validity of the differential pattern (7/8 rounds)

	bool zeroInternalDiffAttack(UINT32 constant[]);//14-round zero-internal-differential distinguisher
private:
	bool constant[6][32];//round constant
	int **threeVarLinearZero;//constant is zero(4*5)
	int **fourVarLinearZero;//constant is zero (8*6)
	int **andConstraint;//9*7
	int **orConstraint;//9*7
	int **andValueConstraint;//12*7
	int **orValueConstraint;//12*7
	char sign[3] = { 'X','Y','Z' };
	UINT64 previous;
	UINT32 previousINT32;
	bool isConstantAddition;
	UINT32 solNum;
};


#endif // !GIMLI_H_
#pragma once
