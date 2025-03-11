//MT-tracker installer
//Bioinformatics Group, Qingdao University
//Updated at Jan. 16, 2025
//Updated by Wenjie Zhu,  Xiaoquan Su

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <sys/stat.h>
#include <unistd.h>
#include <omp.h>
#include <string.h>
#include <numeric>
#include <atomic>
#include <chrono>
#include <thread>

#include "MT_tracker_abd.h"
#include "MT_tracker.h"

using namespace std;

char Ref_db;

string Listfilename;
string Listprefix;
string Queryfile1;
string Queryfile2;
string Virfile;

string Tablefilename;
string Outfilename;
string Outfilename1;
string Outfilename2;
string Ref_name;

_Mt_Tracker_Abd mt_tracker;
int Coren = 0;//Number of threads

bool Is_cp_correct; //Used to specify whether to perform copy correction
bool Is_sim; //Output format distance (T) or similarity (F)
//bool Is_heatmap;
//int Cluster = 2;
bool Is_detail;

int Mode = 0; //0: single, 1: multi_list, 2: multi_table
bool Reversed_table = false;


int printhelp(){
    
    cout << "MT-tracker version : " << Version << endl;
    cout << "\tCompute the transition among samples" << endl;
    cout << "Usage: " << endl;
    cout << "MT-tracker [Option] Value" << endl;
    cout << "Options: " << endl;
    cout << "\t-D (upper) Reference database, default is G (Greangenes13_97%)"<< endl;
    cout << "m (MetaPhlAn2)"<< endl;

    cout << "\t  -i Two samples path for single sample transition" << endl;
    cout << "\tor" << endl;
    cout << "\t  -l Input files list for multi-sample transition" << endl;
    cout << "\t  -p List files path prefix [Optional for -l]" << endl;
    cout << "\tor" << endl;

    cout << "\t  -T (upper) Input species table for multi-sample transition" << endl;
    cout << "\t  -R (upper) If the input table is reversed, T(rue) or F(alse), default is false [Optional for -T]" << endl;
    cout << "\t  -g Input the transition matrix and meta information to obtain the transition probability between groups" << endl;
    cout << "\t[Output options]" << endl;
    cout << "\t  -o Output file" << endl;
    
    cout << "\t[Other options]" << endl;
    cout << "\t  -t Number of threads, default is auto" << endl;
    cout << "\t  -h Help" << endl;
    
    exit(0);
    
    return 0;
    
    };
    
int Parse_Para(int argc, char * argv[]){
    
    Ref_db = 'G';
    Coren = 0;
    Mode = 0; //default is single;
    
    //Is_cp_correct = true;
    Is_sim = true;
    Reversed_table = false;

    int i = 1;
      if (argc ==1) 
		printhelp();
//    If the argument does not begin with a '-', an error message is output and the program exits.
      while(i<argc){
         if (argv[i][0] != '-') {
                           cerr << "Argument # " << i << " Error : Arguments must start with -" << endl;
                           exit(0);
                           };           
         switch(argv[i][1]){
                            case 'D': Ref_db = argv[i+1][0]; break;
                            case 'i': Queryfile1 = argv[i+1]; Queryfile2 = argv[i+2]; i++; Mode = 3; break;

	                        case 'g': Queryfile1 = argv[i+1]; Queryfile2 = argv[i+2];
                                if ((argv[i+3][0] == 'f') || (argv[i+3][0] == 'F')) Is_detail = false;
                                if ((argv[i+3][0] == 't') || (argv[i+3][0] == 'T')) Is_detail = true; i++; Mode = 4; break;

                            case 'l': Listfilename = argv[i+1]; Mode = 1; break;

                            case 'p': Listprefix = argv[i+1]; break;

                            case 'T': Tablefilename = argv[i+1]; Mode = 2; break;

                            case 'R': if ((argv[i+1][0] == 't') || (argv[i+1][0] == 'T')) Reversed_table = true; break;

                            case 'o': Outfilename = argv[i+1]; break;

                            case 't': Coren = atoi(argv[i+1]); break;

                            case 'h': printhelp(); break;

                            default : cerr << "Error: Unrec argument " << argv[i] << endl; printhelp(); break; 
                            }
         i+=2;
         }
    //Get the maximum number of cores available on the system and store it in the max_core_number variable.
    int max_core_number = sysconf(_SC_NPROCESSORS_CONF);
    //If Coren is less than or equal to 0 or greater than max_core_number, Coren is set to max_core_number.
    if ((Coren <= 0) || (Coren > max_core_number)){
                    //cerr << "Core number must be larger than 0, change to automatic mode" << endl;
                    Coren =  max_core_number;
                    } 
    /*
    if (Cluster <= 0){
                cerr << "Warning: cluster number must be larger than 0, change to default (2)" << endl;
                }
    */
return 0;    
}
//Progress bar function
void display_progress_bar(std::atomic<long> &progress, long total, int bar_width = 50) {
    const char postfix[] = {'|', '/', '=', '\\'};
    while (progress < total) {
        int completed = progress * bar_width / total; // Length of completed section
        std::string bar(completed, '#');
        bar.resize(bar_width, ' '); // Fill the unfinished part with blanks

        // Dynamic suffix effect
        int postfix_index = progress % 4;

        // Print progress bar
        std::cout << "\r[" << bar << "] "
                  << std::setw(6) << std::fixed << std::setprecision(2)
                  << (float)progress * 100 / total << "% "
                  << postfix[postfix_index];
        std::cout.flush();

        std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    }

    std::cout << "\r[" << std::string(bar_width, '#') << "] 100.00% " << std::endl;
}

//Used to write matrix data to a file
// n: the size of the matrix, indicating the number of samples
//sim_matrix: a vector pointer to store the similarity matrix
//is_sim: a Boolean value indicating whether to output a similarity matrix or a distance matrix
//sam_name: a vector to store the sample name
void Output_Matrix(const char * outfilename, int n, vector <float> * sim_matrix, bool is_sim, vector <string> sam_name){
    FILE * outfile = fopen(outfilename, "w");
    if (outfile == NULL){
        cerr << "Error: Cannot open output file : " << outfilename << endl;
        return;
    }

    //Label
    for(int i = 0; i < n; i ++)
        fprintf(outfile, "\t%s", sam_name[i].c_str());
    fprintf(outfile, "\n");

    for(int i = 0; i < n; i ++){
        fprintf(outfile, "%s", sam_name[i].c_str());

        for (int j = 0; j < n; j ++){

            long ii = (i <= j) ? i : j;
            long jj = (i >= j) ? i : j;
            long p = ii * (long) n + jj - (1 + ii + 1) * (ii + 1) / 2;

            if (is_sim){
                if (ii == jj|| i > j) fprintf(outfile, "\t0.00");
                else fprintf(outfile, "\t%f", (*sim_matrix)[p]);
            }
            else {
                if (ii == jj||i > j) fprintf(outfile, "\t0.00");
                else fprintf(outfile, "\t%f", 1-(*sim_matrix)[p]);
            }
        }
        fprintf(outfile, "\n");
    }

    fclose(outfile);
}


void Single_Comp2(){
    if(Ref_name.size() == 0)
        mt_tracker = _Mt_Tracker_Abd(Ref_db);

    float * Abd1 = new float [mt_tracker.Get_LeafN()];//Count (number of lines in the file), Abd represents abundance
    float * Abd2 = new float [mt_tracker.Get_LeafN()];//This is a float pointer that will point to the first element of the allocated float array.

    //The function loads the data from Queryfile1 and Queryfile2 and stores them in the corresponding arrays Abd1 and Abd2.
    int mapped_otu_count1 = mt_tracker.Load_abd(Queryfile1.c_str(), Abd1, Is_cp_correct);
    int mapped_otu_count2 = mt_tracker.Load_abd(Queryfile2.c_str(), Abd2, Is_cp_correct);
    if ((mapped_otu_count1) == 0 && (mapped_otu_count2) == 0) {
        cerr << "Error: Failed to map OTUs/Sp for file. Maybe check if database selection is correct." << endl;
        delete[] Abd1;
        delete[] Abd2;
        return;
    }
    tuple<vector<float>, vector<float>, vector<pair<int, int>>> result = mt_tracker.myFunction(Abd1, Abd2);

    vector<float> leaf_abd;//Leaf node abundance
    vector<float> ans_abd;//Ancestral node abundance
    vector<pair<int,int>> vir_id;//Virtual code

    tie(leaf_abd, ans_abd, vir_id)=result;
    std::vector<string> id = mt_tracker.getVector();
    int k = mt_tracker.Get_OredrN();//Number of nodes
    float sum = 0;
    for(int i=0; i< k; i++){
        sum += (i < leaf_abd.size()) ? leaf_abd[i] : 0;
        sum += (i < ans_abd.size()) ? ans_abd[i] : 0;
    }
    sum /= 100.0;
    for(int i=0; i< k; i++){
        if (i < leaf_abd.size()) {
            leaf_abd[i] /= sum;

        }
        if (i < ans_abd.size()) {
            ans_abd[i] /= sum;
        }
    }
    cout <<"virtual ancestor abundance: "<< endl;
    cout << left << setw(50) << "ID" << setw(10) << "Abundance" << endl;
    cout << string(60, '-') << endl;
    //Leaf node abundance and id
    int z = leaf_abd.size();
    for(int i=0; i< z; ++i) {
        if(leaf_abd[i] != 0){
            cout << left << setw(50) << id[i] << setw(10) << leaf_abd[i] << endl;
        }
    }

    cout << endl;

    Is_sim = true;
    vector<vector<float>> ANS_info;//Leaf node abundance and id
    int count = mt_tracker.ANS_count(vir_id, ans_abd, k, ANS_info);//Number of common ancestors
    cout <<"They have " <<count<< " ancestor" <<endl;
    if(count!=0) {
        //Calculate the similarity between Abd1 and ANS
        float sim1 = mt_tracker.Ans_sim(Abd1, leaf_abd, count, ANS_info);
        //Calculate the similarity between Abd2 and ANS
        float sim2 = mt_tracker.Ans_sim(Abd2, leaf_abd, count, ANS_info);
        if (Is_sim)
            cout <<"distance between sample1 and virtual ancestor is "<<1- sim1<<endl;
            cout <<"distance between sample2 and virtual ancestor is "<<1- sim2<<endl;
            if (sim1 < sim2) {
                cout << "sample1 <- sample2" << endl;
            } else if (sim1 > sim2) {
                cout << "sample1 -> sample2" << endl;
            } else {
                cout << "sample1 <-> sample2" << endl; // 如果相等，输出两者平等关系
            }

    }
    if (count==0){
        cerr << "They don't have ancestor "<< endl;
        return;
    }

}

void Multi_Comp(){

     if(Ref_name.size() == 0)
       mt_tracker = _Mt_Tracker_Abd(Ref_db);

     //load list
    vector <string> sam_name;
    vector <string> file_list;
    int file_count = 0;
    file_count = Load_List(Listfilename.c_str(), file_list, sam_name, Listprefix);
            
    //load abd
    const int leafN = mt_tracker.Get_LeafN();
    vector<float> Abd_flat(file_count * leafN);
    vector<float*> Abd(file_count);
    int mapped_otu_count = 0;
    for (int i = 0; i < file_count; i ++){
        Abd[i] = &Abd_flat[i * leafN];
        mapped_otu_count = mt_tracker.Load_abd(file_list[i].c_str(), Abd[i], Is_cp_correct);
        if (mapped_otu_count == 0) {
            cerr << "Error: Failed to map OTUs/Sp for file: " << file_list[i].c_str() << endl;
            cerr << "Maybe check if database selection is correct." << endl;
            // Continue to the next file
            continue;
        }
    }
    if (file_count == 0){
        return;
    }
    cout << file_count << " samples loaded" << endl;
    
    //make order
    vector <int> order_m;
    vector <int> order_n;
    long iter = 0;
    for (int i = 0; i < file_count - 1; i ++)
        for (int j = i + 1; j < file_count; j ++){            
            order_m.push_back(i);
            order_n.push_back(j);
            iter ++;
            }
    //Create and initialize a storage file similarity matrix.
    vector <float>  sim_matrix;
    for (long i = 0; i < iter; i ++)
        sim_matrix.push_back(0);
    // Tracking progress using atomic variables
    std::atomic<long> progress(0);
    // Start the progress bar thread
    std::thread progress_printer([&]() {
        display_progress_bar(progress, iter);
    });
    //openmp
    omp_set_num_threads(Coren);
    cout <<"Using"<< Coren <<"Coren"<< endl;
    int k = mt_tracker.Get_OredrN();
    #pragma omp parallel for schedule(dynamic, 1)
    for (long i = 0; i < iter; i ++){

        long m = order_m[i];
        long n = order_n[i];
        long p = m * (long) file_count + n - (1 + m + 1) * (m + 1) / 2;

        vector<float> leaf_abd;//Leaf node abundance
        vector<float> ans_abd;//Ancestral node abundance
        vector<pair<int,int>> vir_id;//Virtual code

        tie(leaf_abd, ans_abd, vir_id)= mt_tracker.myFunction(Abd[m], Abd[n]);
        float sum = 0;
        for(int i=0; i< k; i++){
            sum += (i < leaf_abd.size()) ? leaf_abd[i] : 0;
            sum += (i < ans_abd.size()) ? ans_abd[i] : 0;
        }
        if (sum != 0){
            sum /= 100.0;
        } else
            sum = 1;
        for(int i=0; i< k; i++){
            if (i < leaf_abd.size()) {
                leaf_abd[i] /= sum;

            }
            if (i < ans_abd.size()) {
                ans_abd[i] /= sum;
            }
        }


        vector<vector<float>> ANS_info;//Public Ancestral node abundance and its child node IDs
        int count = mt_tracker.ANS_count(vir_id, ans_abd, k, ANS_info);//Number of common ancestors

        if(count!=0) {

            //Calculate the similarity between Abd1 and ANS
            float sim1 = mt_tracker.Ans_sim(Abd[m], leaf_abd, count, ANS_info);
            //Calculate the similarity between Abd2 and ANS
            float sim2 = mt_tracker.Ans_sim(Abd[n], leaf_abd, count, ANS_info);
            if (Is_sim)
                sim_matrix[p] = sim1 * sim1 - sim2 * sim2;//Multiply by a weight
            else
                sim_matrix[p] = sim1 - sim2 + 1;
        }

        if(count==0) {

            sim_matrix[p] = 1.0;
        }
        progress++;

    }
    progress_printer.join();
    cout << endl;
    Output_Matrix(Outfilename.c_str(), file_count, &sim_matrix, Is_sim, sam_name);

//    for (int i = 0; i < file_count; i ++)
//        delete [] Abd[i];
//    delete [] Abd;
    

     }

void Multi_Comp_Table(_Table_Format abd_table){

     if(Ref_name.size() == 0)
       mt_tracker = _Mt_Tracker_Abd(Ref_db);

    int file_count = abd_table.Get_Sample_Size();
         
    //load abd
    float **Abd = new float * [file_count];
    int mapped_otu_count = 0;
    for (int i = 0; i < file_count; i ++){
        Abd[i] = new float [mt_tracker.Get_LeafN()];
        //cout << mt_tracker.Load_abd(&abd_table, Abd[i], i, Is_cp_correct) << " OTUs in file " << i + 1 << endl;
        mapped_otu_count = mt_tracker.Load_abd(&abd_table, Abd[i], i, Is_cp_correct);
        if (mapped_otu_count == 0) {
            cerr << "Error: Failed to map OTUs/Sp for file"<< Tablefilename.c_str() << endl;
            cerr << "Maybe check if database selection is correct." << endl;
            return;
        }
    }
    if (file_count == 0){
        return;
    }
    cout << file_count << " samples loaded" << endl;
    
    //make order
    vector <int> order_m;
    vector <int> order_n;
    long iter = 0;
    for (int i = 0; i < file_count - 1; i ++)
        for (int j = i + 1; j < file_count; j ++){            
            order_m.push_back(i);
            order_n.push_back(j);
            iter ++;
            }
        
    vector <float>  sim_matrix;
    for (long i = 0; i < iter; i ++)
        sim_matrix.push_back(0);

    //progress
    std::atomic<long> progress(0);
    std::thread progress_printer([&]() {
        display_progress_bar(progress, iter);
    });
    //openmp

    omp_set_num_threads(Coren);
    int k = mt_tracker.Get_OredrN();
    #pragma omp parallel for schedule(dynamic, 1)
    for (long i = 0; i < iter; i ++){
        
        long m = order_m[i];
        long n = order_n[i];
        long p = m * (long) file_count + n - (1 + m + 1) * (m + 1) / 2;


        vector<float> leaf_abd;//Leaf node abundance
        vector<float> ans_abd;//Ancestral node abundance
        vector<pair<int,int>> vir_id;//Virtual code
        tie(leaf_abd, ans_abd, vir_id) = mt_tracker.myFunction(Abd[m], Abd[n]);

        float sum = 0;
        for(int i=0; i< k; i++){
            sum += (i < leaf_abd.size()) ? leaf_abd[i] : 0;
            sum += (i < ans_abd.size()) ? ans_abd[i] : 0;
        }
        if (sum != 0){
            sum /= 100.0;
        } else
            sum = 1;
        for(int i=0; i< k; i++){
            if (i < leaf_abd.size()) {
                leaf_abd[i] /= sum;
            }
            if (i < ans_abd.size()) {
                ans_abd[i] /= sum;
            }
        }

        vector<vector<float>> ANS_info;//Ancestral node abundance
        int count = mt_tracker.ANS_count(vir_id, ans_abd, k, ANS_info);//Number of common ancestors

        if(count!=0) {
            //Calculate the similarity between Abd1 and ANS
            float sim1 = mt_tracker.Ans_sim(Abd[m], leaf_abd, count, ANS_info);
            //Calculate the similarity between Abd2 and ANS
            float sim2 = mt_tracker.Ans_sim(Abd[n], leaf_abd, count, ANS_info);
            if (Is_sim)
                sim_matrix[p] = sim1 * sim1 - sim2 * sim2;
            else
                sim_matrix[p] = sim1 - sim2 + 1;
        }

        if(count==0) {
            cerr << "They don't have ans " << endl;
        }
        progress++;
        }
    progress_printer.join();
    cout << endl;

    Output_Matrix(Outfilename.c_str(), file_count, &sim_matrix, Is_sim, abd_table.Get_Sample_Names());

    for (int i = 0; i < file_count; i ++)
        delete [] Abd[i];

     };

void Output_direction(const char * outfilename, string& direction){

    std::ofstream outFile(outfilename);
    if (!outFile.is_open()) {
        std::cerr << "Error opening output file." << std::endl;
        return;
    }

    outFile << direction;
    outFile.close();
}

void Group_direction() {
    _Mt_Tracker tracker;
    map<string, vector<string>> meta = tracker.Load_meta(Queryfile2.c_str());
    unordered_map<string, unordered_map<string, float>> table = tracker.Load_ANS(Queryfile1.c_str());

    string output;
    vector<float> all_data;

    for (auto it1 = meta.begin(); it1 != meta.end(); ++it1) {
        const vector<string>& sample1 = it1->second;
        for (auto it2 = meta.begin(); it2 != meta.end(); ++it2) {
            const vector<string>& sample2 = it2->second;
            // Skip identical samples
            if (it1->first == it2->first) {
                continue;
            }
            float data = tracker.calculate_data(sample1, sample2, table);
            all_data.push_back(data);
            output.append("Total for " + it1->first + " and " + it2->first + " " + to_string(data) + "\n");
        }
    }

    // Determine whether normalization is needed
    if (meta.size() > 2) {
        tracker.normalize_and_transform(all_data);

        output.clear();
        size_t index = 0;

        // Build output for normalized data while skipping identical samples
        for (auto it1 = meta.begin(); it1 != meta.end(); ++it1) {
            for (auto it2 = meta.begin(); it2 != meta.end(); ++it2) {
                // Skip identical samples
                if (it1->first == it2->first) {
                    continue;
                }
                output.append("Total for " + it1->first + " and " + it2->first + " " + to_string(all_data[index]) + "\n");
                index++;
            }
        }
    }

    Output_direction(Outfilename.c_str(), output);
}

int main(int argc, char * argv[]){
    
    Parse_Para(argc, argv);         
                  
    switch (Mode){

           case 1: Multi_Comp(); break;
           case 2:{
                   _Table_Format table(Tablefilename.c_str(), Reversed_table);//table为丰度的二维矩阵
                   Multi_Comp_Table(table);
                   break;
                   }
           case 3: Single_Comp2(); break;

           case 4: Group_direction();break;

           default: break;
           }
    return 0;
    }
