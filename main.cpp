#include "Gimli.h"
#include <iostream>
#include <string>
#include <ctime>
using namespace std;

void testCorrectnessOfDiffTrail() {
	UINT32 diffLast[5][12] = {//Invalid part of the 12-round differential trail 
		{0,0,0,0,0,0,0,0,0,0x01008000,0x00000200,0x01000000},
		{0,0,0,0,0,0,0x00010002,0,0,0,0x01040002,0x03008000},
		{0,0,0,0,0,0,0,0x02000400,0x00010002,0x020A0480,0x0A000402,0x0A010000},
		{0x02020104,0,0,0x02000100,0,0,0,0x00080004,0x00020004,0,0x14010430,0x1E081480},
		{0,0x04020804,0x02020104,0,0x00020004,0x02000100,0x00000A00,0x10001800,0x00040008,0xB00A0910,0x02186078,0x3C102900}
	};

	UINT32 diffFirst[5][12] = {//valid part of the 12-round differential trail
		{0x04010100,0,0x02008080,0x80010380,0x40010180,0x40010180,0x06010100,0x02000000,0x03018080,0x80100C00,0x40100400,0x40104400},
		{0,0,0,0x80020080,0x00060080,0x00070480,0,0,0,0x80210180,0x40200080,0x00318400},
		{0,0,0,0x00003100,0x00000100,0x80000980,0,0,0,0x80401180,0x80000180,0x80000980},
		{0,0,0,0,0,0,0,0,0,0x80800100,0x80400000,0x80400080},
		{0,0,0,0,0,0,0,0,0,0x80000000,0x80000000,0x80000000}
	};

	UINT32 diffZong[4][12] = {//The invalid 6-round differential trail
		{0,0,0,0xff898081,0,0,0,0,0,0xff898081,0,0},
		{0,0,0,0,0x80618880,0x81ff8980,0,0,0,0,0x80618880,0x81ff8980 },
		{0,0,0,0x42668080,0xc0400000, 0x00011100,0,0,0,0x42668080,0xc0400000, 0x00011100},
		{0,0,0,0x80010080,0x00402000,0x80400080,0,0,0,0x80010080,0x00402000,0x80400080}
	};

	UINT32 diff_6R[5][12] = {//Our found valid 6-round differential trail
		{0,0,0,0x7c2c642a,0,0,0,0,0,0x7c2c642a,0,0},
		{0,0,0,0,0x6e1c342c,0x2a7c2c64,0,0,0,0,0x6e1c342c,0x2a7c2c64},
		{0,0,0,0x91143078,0x28785014,0x35288a58,0,0,0,0x91143078,0x28785014,0x35288a58},
		{0,0,0,0x80010008,0x2000,0x44400080,0,0,0,0x80010008,0x2000,0x44400080},
		{0,0,0,0x80,0x400000,0x80000000,0,0,0,0x80,0x400000,0x80000000}
	};

	Gimli g;
	cout << "Our found 6-round differential trail:" << endl;
	g.testCorrectOfDiff(diff_6R, 5, 24, 19);
	cout << endl << endl;

	cout << "The 6-round differential trail found by Zong et al.:" << endl;
	g.testCorrectOfDiff(diffZong, 4, 24, 20);
	cout << endl << endl;

	UINT32 diff[5][12];
	for (int i = 0; i < 5; i++) {
		for (int j = 0; j < 12; j++) {
			diff[i][j] = diffLast[i][j];
		}
	}
	for (int r = 0; r < 5; r++) {
		for (int j = 0; j < 12; j++) {
			if (j % 3 == 1) {
				diff[r][j] = g.RR(diff[r][j], 9);
			}
			else if (j % 3 == 0) {
				diff[r][j] = g.RR(diff[r][j], 24);
			}
		}
	}
	cout << "The invalid part of the 12-round differential trail:" << endl;
	g.testCorrectOfDiff(diff, 4, 16, 12);
	cout << endl << endl;

	for (int i = 0; i < 5; i++) {
		for (int j = 0; j < 12; j++) {
			diff[i][j] = diffFirst[i][j];
		}
	}
	for (int r = 0; r < 5; r++) {
		for (int j = 0; j < 12; j++) {
			if (j % 3 == 1) {
				diff[r][j] = g.RR(diff[r][j], 9);
			}
			else if (j % 3 == 0) {
				diff[r][j] = g.RR(diff[r][j], 24);
			}
		}
	}
	cout << "The valid part of the 12-round differential trail:" << endl;
	g.testCorrectOfDiff(diff, 5, 24, 19);
}

void verifySolForS1() {
	Gimli g;
	g.verifySolForS1();
}

void zero_internal_differential_test() {
	Gimli g;
	UINT32 constant[4] = { 0xfed8903,0x3456789,0xffdd223,0xf2343dde };
	UINT32 testTimes = EXP[20];
	UINT32 succeedTimes = 0;
	cout << hex << "0x"<<testTimes << " times of experiments..." << endl;
	for (UINT32 i = 0; i < testTimes; i++) {
		constant[0] = g.getRand32();
		constant[1] = g.getRand32();
		constant[2] = g.getRand32();
		constant[3] = g.getRand32();
		if (g.zeroInternalDiffAttack(constant)) {
			succeedTimes++;
		}
	}
	cout << hex << "Testing times:" << hex << testTimes << endl;
	cout << hex << "Successful times:" << hex << succeedTimes << endl;
}

void invalid7RoundDiffPattern() {
	Gimli g;
	g.searchCollision3RPattern();
}

void check8RoundDiffPattern() {
	Gimli g;
	g.searchCollision4RPattern();
}

void searchValid6RoundDiff() {
	Gimli g;
	cout << "Please input the number of threads:";
	int threadNum;
	cin >> threadNum;
	g.search6RSFSCollision(threadNum);
}

void generateValidCapacityPartOfS0() {
	Gimli g;
	g.generateValidCapacityFor6RCollision();
	cout << endl<<"Count the number of unique values" << endl;
	g.countValidCapacityFor6RCollision();
}

void verifyInstances() {
	Gimli g;
	g.checkRealSFSCollision();
}

int main(){
	cout << "0 -> Check the correctness of the diff trail (valid 6-round trail, invalid 6-round trail, 12-round trail)" << endl;
	cout << "1 -> Check the correctness of the solutions of S^1 in our 6-round collision attack" << endl;
	cout << "2 -> Test the 14-round zero-internal-differential distinguisher" << endl;
	cout << "3 -> The invalid differential pattern for 7-round SFS collision attack" << endl;
	cout << "4 -> Gurobi-solver-based verification of our diff pattern for 8-round SFS collision attack" << endl;
	cout << "5 -> Search the 6-round diff trail and conforming message pair simultaneously" << endl;
	cout << "6 -> Compute the valid solutions of the capacity part of S^0 in our collision attack" << endl;
	cout << "7 -> Verify the SFS colliding message pair in the paper" << endl;
	cout << endl<<"Please input your command:";
	
	int cmd;
	cin >> cmd;
	if (cmd == 0) {
		testCorrectnessOfDiffTrail();
	}
	else if (cmd == 1) {
		verifySolForS1();
	}
	else if (cmd == 2) {
		zero_internal_differential_test();
	}
	else if (cmd == 3) {
		invalid7RoundDiffPattern();
	}
	else if (cmd == 4) {
		check8RoundDiffPattern();
	}
	else if (cmd == 5) {
		searchValid6RoundDiff();
	}
	else if (cmd == 6) {
		generateValidCapacityPartOfS0();
	}
	else if (cmd == 7) {
		verifyInstances();
	}
	else {
		cout << "wrong command!" << endl;
	}
	return 0;
}