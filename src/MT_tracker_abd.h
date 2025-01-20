//MT-tracker installer
//Bioinformatics Group, Qingdao University
//Updated at Jan. 16, 2025
//Updated by Wenjie Zhu,  Xiaoquan Su

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <vector>

#include "hash.h"
#include "utility.h"
#include "version.h"
#include "table_format.h"
#include "db.h"
#include "otu_parser.h"
#include "sp_parser.h"
#include "MT_tracker.h"
#include "dist.h"

#ifndef MT_TRACKER_ABD_H
#define MT_TRACKER_ABD_H

using namespace std;

class _Mt_Tracker_Abd : public _Mt_Tracker {
      
      public:
             _Mt_Tracker_Abd() : _Mt_Tracker() {};
    
            _Mt_Tracker_Abd(char db) : _Mt_Tracker(db){
                        
                        //Database.Set_DB(db);                                                  
                        Sp_parser = _SP_Parser(Database);                                             
                        }

    
             int Load_abd(const char * infilename, float * Abd, bool is_cp_correct);
             int Load_abd(const char * infilename, float * Abd);
             int Load_ANS_abd(const char * infilename, float * Abd, bool is_cp_correct);
             int Load_abd(_Table_Format * table, float * Abd, int sample, bool is_cp_correct); //Load by table_format
             int Load_abd(_Table_Format * table, float * Abd, int sample); //Load by table_format


      private:
               _SP_Parser Sp_parser;
    
              };
void ConvertHashMap(const hash_map<string, int, std_string_hash> &src, hash_map<string, float, std_string_hash> &dest) {
    for (auto &entry : src) {
        dest[entry.first] = static_cast<float>(entry.second);
    }
}
int _Mt_Tracker_Abd::Load_ANS_abd(const char * infilename, float * Abd, bool is_cp_correct){

    memset(Abd, 0, LeafN * sizeof(float));
    /*These three lines of code define three hash maps, which are used to store data for taxon_count, otu_count, and otu_abd.
     * These hash maps use strings as keys and floating point numbers as values.
temp_taxon_count is to correspond to the function in the pms file*/
    hash_map<string, int, std_string_hash> temp_taxon_count;
    hash_map<string, float, std_string_hash> taxon_count;
    Otu_parser.Load_file_to_hash(infilename, temp_taxon_count);
    ConvertHashMap(temp_taxon_count, taxon_count);
    hash_map<string, float, std_string_hash> otu_count;
    hash_map<string, float, std_string_hash> otu_abd;
    //Calculate and populate the otu_count hash map based on the data in the taxon_count hash map.
    Sp_parser.Load_sp_to_hash(taxon_count, otu_count);


    //cp_number_correct
    for (hash_map<string, float, std_string_hash>::iterator miter = otu_count.begin(); miter != otu_count.end(); miter ++){

        float cp_no = 1.0;
        if (is_cp_correct)
            cp_no = Otu_parser.Get_cp_by_OTU(miter->first);

        otu_abd[miter->first] = (float) miter->second / cp_no;

    }

    //Initialize a variable mapped_otu_count to 0 to calculate the number of OTUs mapped to the Abd array.
    int mapped_otu_count = 0;
    for (hash_map<string, float, std_string_hash>::iterator miter = otu_abd.begin(); miter != otu_abd.end(); miter ++){

        //cout << miter->first << endl;
        //If the OTU exists in the Id_hash mapping, its relative abundance is stored in the corresponding position of the Abd array and the mapped_otu_count count is increased.
        if (Id_hash.count(miter->first) != 0){
            Abd[Id_hash[miter->first]] = miter->second;
            mapped_otu_count ++;
        }
    }
    //Returns mapped_otu_count, which indicates the number of OTUs successfully mapped to the Abd array.
    return mapped_otu_count;
}


int _Mt_Tracker_Abd::Load_abd(const char * infilename, float * Abd, bool is_cp_correct){

    memset(Abd, 0, LeafN * sizeof(float));
    hash_map<string, int, std_string_hash> temp_taxon_count;
    hash_map<string, float, std_string_hash> taxon_count;
    Otu_parser.Load_file_to_hash(infilename, temp_taxon_count);
    ConvertHashMap(temp_taxon_count, taxon_count);
    hash_map<string, float, std_string_hash> otu_count;
    hash_map<string, float, std_string_hash> otu_abd;
    Sp_parser.Load_sp_to_hash(taxon_count, otu_count);
    
    float total = 0;
        
    //cp_number_correct
    for (hash_map<string, float, std_string_hash>::iterator miter = otu_count.begin(); miter != otu_count.end(); miter ++){
        
        float cp_no = 1.0;
        if (is_cp_correct)
           cp_no = Otu_parser.Get_cp_by_OTU(miter->first);
           
        otu_abd[miter->first] = (float) miter->second / cp_no;
        total += otu_abd[miter->first];
        }
        
    //norm
    total /= 100.0;
    int mapped_otu_count = 0;
    for (hash_map<string, float, std_string_hash>::iterator miter = otu_abd.begin(); miter != otu_abd.end(); miter ++){
        miter->second /= total;
        //debug
        //cout << miter->first << endl;
        if (Id_hash.count(miter->first) != 0){
            Abd[Id_hash[miter->first]] = miter->second;
            mapped_otu_count ++;
            }
        }
    return mapped_otu_count;
    }

int _Mt_Tracker_Abd::Load_abd(const char * infilename, float * Abd){
    
    return this->Load_abd(infilename, Abd, true);
    }

int _Mt_Tracker_Abd::Load_abd(_Table_Format * table, float *Abd, int sample, bool is_cp_correct){
    
    memset(Abd, 0, LeafN * sizeof(float));

    vector <string> otus = table->Get_Feature_Names();
    vector <float> abds = table->Get_Abd(sample);
    //cout<< otus[1] <<" "<<abds[1]<<endl;
    hash_map<string, float, std_string_hash> taxon_count;
    hash_map<string, float, std_string_hash> otu_count;
    hash_map<string, float, std_string_hash> otu_abd;
    
    //load into hash
    for (int i = 0; i < otus.size(); i ++)
        if (abds[i] > 0){
        //Check whether the naming of input samples is standardized
           //string a_otu = Check_SP(otus[i]);
           string a_otu = otus[i];
           if (taxon_count.count(a_otu) == 0) taxon_count[a_otu] = abds[i];
           else taxon_count[a_otu] += abds[i];
           }
        
    Sp_parser.Load_sp_to_hash(taxon_count, otu_count); 
    
    float total = 0;
    
    //cp_number_correct
    for (hash_map<string, float, std_string_hash>::iterator miter = otu_count.begin(); miter != otu_count.end(); miter ++){
        
        float cp_no = 1.0;        
        if (is_cp_correct)
           cp_no = Otu_parser.Get_cp_by_OTU(miter->first);
           
        otu_abd[miter->first] = (float) miter->second / cp_no;
        total += otu_abd[miter->first];
        }
    //norm
    total /= 100.0;
    int mapped_otu_count = 0;
    for (hash_map<string, float, std_string_hash>::iterator miter = otu_abd.begin(); miter != otu_abd.end(); miter ++){
                          miter->second /= total;
			  //cout<<"otu:"<< miter->first<<",abd:"<<miter->second<<endl;
                          if (Id_hash.count(miter->first) != 0){
                             Abd[Id_hash[miter->first]] = miter->second;
                             mapped_otu_count ++;
                             }
                          }
    //cout<<mapped_otu_count<<endl;
    return mapped_otu_count;
    }

int _Mt_Tracker_Abd::Load_abd(_Table_Format * table, float *Abd, int sample){
                                       return this->Load_abd(table, Abd, sample, true);
                                       }

#endif
