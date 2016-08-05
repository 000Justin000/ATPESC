/**
 *	
 *	Charm Tutorial Exercise
 *	key-value store, reduction, sdag
 *
 *	
 **/

#include <vector>
#include <map> 
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "charm++.h"
#include "main.decl.h"

using namespace std;
/*readonly*/ CProxy_Main mainProxy;
/*readonly*/ CProxy_KeyValueStore kvstoreProxy;
/*readonly*/ CProxy_KeyValueClient kvclientProxy;
/*readonly*/ int N; //number of elements in KeyValueStore array
/*readonly*/ int M; //keys stored per array element
/*readonly*/ int K; //number of key requests from an array element

struct KeyValue
{
	int key;
	int value;
};

class Main : public CBase_Main
{
public:
    Main(CkArgMsg *m)
    {
        mainProxy = thisProxy;
        N = CkNumPes();
        M = 10;
        K = 3;
			
		ckout << "---------------------------------------------------------------------" << endl;
		ckout << "Num Array Elements (N): "                                         << N << endl;
		ckout << "Num Keys Stored Per Array Element of KeyValueStore (M): "         << M << endl;
		ckout << "Num Keys Requested By Each Array Element of KeyValueClient (K): " << K << endl;
		ckout << "caling method send_request in KeyValueClient array ...."               << endl;
		ckout << "---------------------------------------------------------------------" << endl;
        kvstoreProxy = CProxy_KeyValueStore::ckNew(N);
		kvclientProxy = CProxy_KeyValueClient::ckNew(N);
		kvclientProxy.send_requests();
    };

    void finish()
    {
	  	ckout << "all responses received ... exiting ..." << endl;
        CkExit();	
    }
};



class KeyValueStore : public CBase_KeyValueStore
{
private:
    std::map<int, int> kvmap;

public:
    KeyValueStore()
    {
      for(int i=0; i<M; i++)
        kvmap[thisIndex*M + i] = thisIndex*M + i;
    }

    KeyValueStore(CkMigrateMessage* m) {};

	void request(int refnum, int key, int reqIdx)
    {
        printf("requesting from client[%d] to store [%d]\n", reqIdx, thisIndex);
		kvclientProxy[reqIdx].response(refnum, kvmap[key], thisIndex);
	}
};



class KeyValueClient : public CBase_KeyValueClient
{
KeyValueClient_SDAG_CODE
private:
    KeyValue* kvpairs;
    int num_finished;

public:
    KeyValueClient()
    {
        kvpairs = new KeyValue[K];
    };

    KeyValueClient(CkMigrateMessage* m) {};

    void response(int refnum, int value, int resIdx)
    {
        printf("responsing from store [%d] to client[%d]: %d\n", resIdx, thisIndex, value);
        kvpairs[refnum].value = value;      
        num_finished++; 

        if (num_finished == K)
        {
            printf("client[%d]: job finished\n", thisIndex);
            contribute(CkCallback(CkReductionTarget(Main, finish), mainProxy)); 
        }
    }

    void send_requests()
    {
	    for(int i=0; i<K; i++)
        {
	    	kvpairs[i].key = i*M;
	    }
        num_finished = 0;

        for(int i=0; i<K; i++)
        {
		    kvstoreProxy[i].request(i, kvpairs[i].key, thisIndex);
        }

        printf("client[%d]: waiting for response\n", thisIndex);
	};
};

#include "main.def.h"
