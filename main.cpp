/*
 * main.cpp
 *
 *  Created on: 14 Oct 2022
 *      Author: Xinjie ZHOU
 */
#include "head.h"

int main(int argc, char** argv){

    if( argc < 4 || argc > 11){//
        printf("usage:\n<arg1> source path, e.g. /export/project/xzhouby\n");
        printf("<arg2> name of dataset, e.g. NY\n");
        printf("<arg3> tree width, e.g. 20\n");
        printf("<arg4> (optional) algorithm for core construction, (0: BPCL; 1: PCL; 2: PLL; 3: WPSL; 4: GLL; 5: Read) default: 0\n");
        printf("<arg5> (optional) algorithm for core update, (0: PDPLL; 1: SDPLL), default: 0\n");
        printf("<arg6> (optional) update type, (0: No Update Test; 1: Decrease; 2: Increase), default: 0\n");
        printf("<arg7> (optional) thread number, default: 15\n");
        printf("<arg8> (optional) batch size for BPCL, (>= thread number), default: thread number\n");
        printf("<arg9> (optional) tree index strategy, (0: local tree index; 1: global tree index), default: 0\n");
        printf("<arg10> (optional) if use extended label, (0: false; 1: True), default: 0\n");
        exit(0);
    }

    string DesFile="./data/";
    string dataset = "NY";
    int treeWidth = 20;
    int algoCoreC = 0;
    int algoCoreU = 0;
    int updateType = 0;
    int runtimes = 1000;
    int updateBatch = 10;
//    updateBatch = 1000;
    updateBatch = 50;
    int threadNum = 15;
    int batchSize = 15;

    bool ifOpt = false;//if use query-oriented optimization
//    ifOpt = true;
    bool ifExtension = false;//if use extension optimization

    if(argc > 1) {
        cout << "argc: " << argc << endl;
        cout << "argv[1]: " << argv[1] << endl;//source path
        DesFile = argv[1];
        if(argc > 2){
            cout << "argv[2]: " << argv[2] << endl;//dataset
            dataset = argv[2];
        }
        if(argc > 3){
            cout << "argv[3]: " << argv[3] << endl;//treewidth
            treeWidth = stoi(argv[3]);
        }
        if(argc > 4){
            cout << "argv[4]: " << argv[4] << endl;//algorithm for core construction
            algoCoreC = stoi(argv[4]);
        }
        if(argc > 5){
            cout << "argv[5]: " << argv[5] << endl;//algorithm for core update
            algoCoreU = stoi(argv[5]);
        }
        if(argc > 6){
            cout << "argv[6]: " << argv[6] << endl;//update type
            updateType = stoi(argv[6]);
        }
        if(argc > 7){
            cout << "argv[7]: " << argv[7] << endl;//thread number
            threadNum = stoi(argv[7]);
        }
        if(argc > 8){
            cout << "argv[8]: " << argv[8] << endl;//batch size
            batchSize = stoi(argv[8]);
        }
        if(argc > 9){
            cout << "argv[9]: " << argv[9] << endl;//if used query-orient optimization
            ifOpt = stoi(argv[9]);
        }
        if(argc > 10){
            cout << "argv[10]: " << argv[10] << endl;//if used extension optimization
            ifExtension = stoi(argv[10]);
        }
    }


	//used for running time calculation
    Timer tt0;
    tt0.start();

//    string graphfile="/media/TraminerData/mengxuan/MengxuanGraphWPSL/Cond/CondWeighted";
    string graphfile=DesFile+"/"+dataset+"/"+dataset;
//    string orderfile=graphfile+".order";
    string ODfile=graphfile+".query";
    string updateFile=graphfile+".update";

    Graph g;
    g.threadnum=threadNum;//thread number of parallel computation (can be changed)
    if(batchSize>=threadNum){
        g.batchsize=batchSize;
    }else{
        g.batchsize=threadNum;
    }

    g.ifParallel = true;
    g.bandWidth=treeWidth;//15;
    g.dataset=dataset;
    g.ifOpt=ifOpt;
    g.ifExtension=ifExtension;
    g.algoCoreC=algoCoreC;
    g.algoCoreU=algoCoreU;
    cout<<"Dataset: "<<dataset<<endl;
    cout<<"Bandwidth: "<<g.bandWidth<<endl;
    cout<<"Thread number: "<<g.threadnum<<endl;
    cout<<"Batch size: "<<g.batchsize<<endl;

    cout<<"This is test for DCT!"<<endl;
    if(ifOpt){
        cout<<"If query-orient optimization: Yes"<<endl;
    }else{
        cout<<"If query-orient optimization: No"<<endl;
    }
    if(ifExtension){
        cout<<"If extension optimization: Yes"<<endl;
    }else{
        cout<<"If extension optimization: No"<<endl;
    }

    if(g.algoCoreU==0){
        cout<<"Core Update algorithm: Propagation-based DPLL."<<endl;
    }else if(g.algoCoreU==1){
        cout<<"Core Update algorithm: Search-based DPLL."<<endl;
    }


//    g.CoreGraphDebug(graphfile+"Core3");//Test
//    g.CoreGraphDebug(graphfile+"CorePLL_148");//NY, original

    g.graphfile=graphfile;


    g.ReadGraph(graphfile);//
//    g.StainingMethod(0);


    ///Task 1: Index construction
    g.CTIndexConstruct();
//    g.WriteCTIndex(graphfile);



    ///Task 2: Query processing
    g.CorrectnessCheck(100);
    g.EffiCheck(ODfile,runtimes);//query efficiency test
//    g.SameTreeQueryTest(ODfile,runtimes);
//    exit(0);
    ///Task 3: Index update
    g.IndexMaintenance(updateFile+"2",updateType,updateBatch);//index maintenance
//    g.IndexMaintenance(updateFile+"ST",updateType,updateBatch);//same-tree index maintenance

    tt0.stop();
    cout<<"\nOverall runtime: "<<tt0.GetRuntime()<<" s."<<endl;
    cout<<"------------------\n"<<endl;
    exit(0);
    g.clear();
	return 0;
}
