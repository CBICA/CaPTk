#include "LBPFeatures2D.h"

LBPFeatures::LBPFeatures(){

}
LBPFeatures::~LBPFeatures(){

}
unsigned int rotateLeft(unsigned int i, unsigned int samples) {
	unsigned int bg = ((i & (1 << (samples - 1))) >> (samples - 1)); // bitget(r,samples)
	unsigned int bs = (i << 1) & ((int)pow(2., (int)samples) - 1); // bitshift(r, 1, samples)
	unsigned int j = (bs + bg) & ((int)pow(2., (int)samples) - 1); // bitset( bs, 1, bg )
	return j;
}
int NumberOfSetBits(int i) {
	i = i - ((i >> 1) & 0x55555555);
	i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
	return (((i + (i >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}
std::vector<float> LBPFeatures::hist(std::vector<float>& src, int bins)
{

	std::vector<float> hist(bins, 0);
	for (float i = 0; i < src.size(); i++) {
		int bin = src.at(i);
		hist[bin] += 1;
	}
	for (float i = 0; i < hist.size(); i++) {

		hist[i] = hist[i] / src.size();
	//	std::cout << hist[i];
	}


	return hist;
}

LBPFeatures::mapping_table LBPFeatures::LUT(unsigned int samples, int type)
{
	std::vector<int> table[256];
	int newMax = 0; //number of patterns in the resulting LBP code
	int index = 0;
	if (type == 0) {

		newMax = static_cast<int>(pow(2., (int)samples));
		for (int i = 0; i < newMax; i++) {
			table->push_back(i);
		}
	}

	else if (type == 1) {
		// Uniform 2
		newMax = samples * (samples - 1) + 3;

		for (unsigned int i = 0; i < pow(2., (int)(samples)); i++) {
			unsigned int j = rotateLeft(i, samples);

			int numt = NumberOfSetBits(i ^ j);

			if (numt <= 2) {
				table->push_back(index);
				index = index + 1;
			}
			else {
				table->push_back(newMax - 1);
			}
		}
	}
	else if (type == 2) {
		long N = (int)pow(2., (int)samples);
		// Rotation Invariant
		int * tmpMap = new int[N];
		memset((void *)tmpMap, -1, N);

		for (long i = 0; i < N; i++) {
			tmpMap[i] = -1;

			unsigned long rm = i;
			unsigned long r = i;
			for (unsigned int j = 1; j <= samples - 1; j++) {
				r = rotateLeft(r, samples);
				if (r < rm)
					rm = r;
			}
			if (tmpMap[rm] < 0) {
				tmpMap[rm] = newMax;
				newMax = newMax + 1;
			}
			table->push_back(tmpMap[rm]);
	//		std::cout << table->at(rm);
		}

	}
	else if (type == 3) {
		// Rotation invariant uniform 2
		newMax = samples + 2;
		for (unsigned int i = 0; i <= pow(2., (int)samples) - 1; i++) {
			unsigned int j = rotateLeft(i, samples); //bitset( bitshift( i, 1, samples ), 1, bitget( i, samples ) ); // rotate left
			unsigned int numt = NumberOfSetBits(i ^ j); //sum(bitget(bitxor(i,j),1:samples));
			if (numt <= 2){
				table->push_back(NumberOfSetBits(i));
		//		std::cout << table->at(i);
			}
			else
			{
				table->push_back(samples + 1);
			//	std::cout << table->at(i);
			}
		}
	}

	mapping_table map_table;
	map_table.N = newMax;
	map_table.table = *table;
	return map_table;
}


