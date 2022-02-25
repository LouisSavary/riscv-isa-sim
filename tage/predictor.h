#ifndef _PREDICTOR_H_
#define _PREDICTOR_H_

#include <inttypes.h>
#include "utils.h"
// #include "predictor.cc"

#include <vector>

class folded_history
{
public:
	unsigned comp;
	int CLENGTH;
	int OLENGTH;
	int OUTPOINT;

	folded_history();

  void init(int original_length, int compressed_length);
  void update(uint8_t *h, int PT);
};


class bentry // TAGE bimodal table entry
{
public:
	int8_t hyst;
	int8_t pred;

	bentry();
};

class gentry // TAGE global table entry
{
public:
	int8_t ctr;
	uint tag;
	int8_t u;

	gentry();
};

class lentry // loop predictor entry
{
public:
	uint16_t NbIter;			// 10 bits
	uint8_t confid;				// 4bits
	uint16_t CurrentIter; // 10 bits

	uint16_t TAG; // 10 bits
	uint8_t age;	// 4 bits
	bool dir;			// 1 bit

	// 39 bits per entry
	lentry();
};

class PREDICTOR
{
public:
	int THRES;

	PREDICTOR(void);
	bool GetPrediction(UINT64 PC, bool isprepred=false);
	void UpdatePredictor(UINT64 PC, OpType opType, bool resolveDir, bool predDir, UINT64 branchTarget);
	void TrackOtherInst(UINT64 PC, OpType opType, bool taken, UINT64 branchTarget);											 
private:
	void reinit();
	int bindex(UINT64 PC);
  int F(long long A, int size, int bank);
	int gindex(unsigned int PC, int bank, long long hist,
						 folded_history *ch_i);
	uint16_t gtag(unsigned int PC, int bank, folded_history *ch0,
								folded_history *ch1);
  void ctrupdate(int8_t &ctr, bool taken, int nbits);
	bool getbim();
	void baseupdate(bool Taken);
	int MYRANDOM();
	void Tagepred(UINT64 PC, bool isprepred);

	void HistoryUpdate(UINT64 PC, OpType opType, bool taken,
										 UINT64 target, long long &X, int &Y,
										 folded_history *H, folded_history *G,
										 folded_history *J);
  int Gpredict(UINT64 PC, long long BHIST, int *length,
							 int8_t **tab, int NBR, int logs, int8_t *W);
  void Gupdate(UINT64 PC, bool taken, long long BHIST, int *length,
							 int8_t **tab, int NBR, int logs, int8_t *W);
  int lindex(UINT64 PC);
  bool getloop(UINT64 PC);
  void loopupdate(UINT64 PC, bool Taken, bool ALLOC);


};


#endif