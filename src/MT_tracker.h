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

#include <unordered_map>

#include "hash.h"
#include "utility.h"
#include "version.h"
#include "table_format.h"
#include "db.h"
#include "otu_parser.h"
#include "dist.h"

#ifndef MT_TRACKER_H
#define MT_TRACKER_H

#define REG_SIZE 150


using namespace std;

class _Mt_Tracker{
      
      public:
             _Mt_Tracker(){
                           LeafN = 0;
                           OrderN = 0;                 
                           //Init();
                          }
    
            _Mt_Tracker(char db){
                        Database.Set_DB(db);                          
                        LeafN = 0;
                        OrderN = 0;                        
                        Init();   
                        Otu_parser = _OTU_Parser(Database);                   
                        }



             tuple<vector<float>, vector<float>, vector<pair<int, int>>> myFunction(float * Abd_1, float * Abd_2);
             int ANS_count(vector<pair<int, int>> &vir_id, vector<float> &ans_abd, int k, vector<vector<float>> &ANS_info);

             float Ans_sim(float * Abd_1, vector<float> Abd_2, int j, vector<vector<float>> ANS_info);
             vector<string>  getVector();
             std::map<std::string, std::vector<std::string>> Load_meta(const char* filepath);
             unordered_map<string, unordered_map<string, float>> Load_ANS(const char * infilename);
             float calculate_data(const vector<string>& sample1, const vector<string>& sample2,
                         const unordered_map<string, unordered_map<string, float>>& table);
             void normalize_and_transform(std::vector<float>& data);


             int Get_LeafN(){
                 return LeafN;
                 }
             int Get_OredrN(){
                 return OrderN;
             }
             int Get_VirN(){
                return VirN;
             }
                 
             string Get_Id(int n){
                    if ((n < 0) || (n >= LeafN)) return "NULL";
                    return Id[n];
                    }
             
      protected:
              _PMDB Database;
              _OTU_Parser Otu_parser;
    
              vector <float> Dist_1;
              vector <float> Dist_2;
              
              vector <int> Order_1;
              vector <int> Order_2;
              vector <int> Order_d;

              vector <string> Id;

              hash_map <string, int, std_string_hash> Id_hash;

              void Init();
              int Load_id();
              int Load_order();


              int LeafN;
              int OrderN;
              int VirN;


};
//The main purpose of this function is to initialize a tree.
// It loads a set of IDs and builds a hash table to find the index based on the ID.
// At the same time, it also checks whether the order of loading the tree needs to be performed and performs corresponding operations.
void _Mt_Tracker::Init(){

     LeafN = Load_id();

     for (int i = 0; i < LeafN; i ++)
        Id_hash[Id[i]] = i;

     OrderN = 0;//Count (how many bacteria are there in the file)
     //load tree  
     if (Database.Get_Is_Tree())      
        OrderN = Load_order();

     }

vector<string> _Mt_Tracker:: getVector() {

    ifstream infile(Database.Get_Tree_Id().c_str(), ifstream::in);
    if (!infile){
        cerr << "Error: Cannot open file : " << Database.Get_Tree_Id() << endl;
//        return false;
    }

    string buffer;
    int count = 0;
    while(getline(infile, buffer)){
        if (buffer.size() == 0) continue;
        Id.push_back(buffer);
        count ++;
    }

    infile.close();
    infile.clear();
    return Id;
}


// Load the calculated data and store it in the table
unordered_map<string, unordered_map<string, float>> _Mt_Tracker::Load_ANS(const char* filepath) {
    ifstream infile(filepath, ifstream::in);
    if (!infile) {
        cerr << "Error: Open Taxonomy Annotation file error" << endl;
        exit(0);
    }

    unordered_map<string, unordered_map<string, float>> matrix;
    string line;
    vector<string> columnNames;

    // 读取文件的第一行，获取列名称
    if (getline(infile, line)) {
        istringstream iss(line);
        string columnName;
        while (iss >> columnName) {
            columnNames.push_back(columnName);
        }
    }

    // 逐行读取数据，填充嵌套哈希表
    while (getline(infile, line)) {
        istringstream iss(line);
        string rowName;
        iss >> rowName;
        float data;
        for (const string& columnName : columnNames) {
            iss >> data;
            matrix[rowName][columnName] = data;
        }
    }

    infile.close();
    infile.clear();

    return matrix;
}

std::string toLower(std::string str) {
    std::transform(str.begin(), str.end(), str.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    return str;
}

std::map<std::string, std::vector<std::string>> _Mt_Tracker::Load_meta(const char* filepath) {

    std::ifstream file(filepath, std::ifstream::in);
    if (!file) {
        std::cerr << "Error: Open Meta file error" << std::endl;
        exit(0);
    }
    // Create a map to associate the numbers with the corresponding vectors
    std::map<std::string, std::vector<std::string>> dataMap;
    std::string buffer;
    std::string SampleID;
    std::string Group;
    getline(file, buffer); // title

    while (getline(file, buffer)) {
        std::istringstream iss(buffer);
        if (getline(iss, SampleID, '\t') && getline(iss, Group, '\t')) {
            // Store the data in the corresponding vector according to the number (ignore case)
            dataMap[toLower(Group)].push_back(SampleID);
        } else {
            std::cerr << "Load meta error" << std::endl;
        }
    }

    file.close();
    file.clear();

    return dataMap;
}

float _Mt_Tracker:: calculate_data(const vector<string>& sample1, const vector<string>& sample2,
                     const unordered_map<string, unordered_map<string, float>>& table) {
    int count = 0;
    float data = 0;

    if (&sample1 == &sample2) {
        for (auto itSample1 = sample1.begin(); itSample1 != sample1.end(); ++itSample1) {
            for (auto itSample2 = itSample1; itSample2 != sample2.end(); ++itSample2) {
                if (table.find(*itSample1) != table.end() &&
                    table.at(*itSample1).find(*itSample2) != table.at(*itSample1).end()) {
                    data += table.at(*itSample1).at(*itSample2);
                    count++;
                }
            }
        }
    } else {
        for (const string& SampleID1 : sample1) {
            for (const string& SampleID2 : sample2) {
                if (table.find(SampleID1) != table.end() && table.at(SampleID1).find(SampleID2) != table.at(SampleID1).end()) {
                    if (table.at(SampleID1).at(SampleID2) == 0) {
                        if (table.at(SampleID2).find(SampleID1) != table.at(SampleID2).end() && table.at(SampleID2).at(SampleID1) != -1) {
                            data += -table.at(SampleID2).at(SampleID1);
                            count++;
                        }
                    } else if (table.at(SampleID1).at(SampleID2) != -1) {
                        data += table.at(SampleID1).at(SampleID2);
                        count++;
                    }
                }
            }
        }
    }

    return (count > 0) ? (data / count) : 0;
}

void _Mt_Tracker::normalize_and_transform(std::vector<float>& data) {
    // 获取最大值和最小值
    float data_max = *std::max_element(data.begin(), data.end());
    float data_min = *std::min_element(data.begin(), data.end());
    float e = 0.01;

    // 保存原始符号
    std::vector<int> signs(data.size());
    for (size_t i = 0; i < data.size(); ++i) {
        signs[i] = (data[i] >= 0) ? 1 : -1;  // 记录符号
        data[i] = std::abs(data[i]);        // 取绝对值
    }

    // 归一化数据
    for (float& value : data) {
        value = e + (1 - 2 * e) * (value - data_min) / (data_max - data_min);
    }

    // 对数变换
    for (float& value : data) {
        value = 1 - log(value + 1);
    }

    // 恢复符号
    for (size_t i = 0; i < data.size(); ++i) {
        data[i] *= signs[i];
    }
}


int _Mt_Tracker:: ANS_count(vector<pair<int, int>> &vir_id, vector<float> &ans_abd, int k, vector<vector<float>> &ANS_info) {

    int j = 0;
    for (int i = 0; i <= k+1; i++) {
        if (ans_abd[i] != 0 && vir_id[i].first < OrderN) {
	    //cout<<ans_abd[i]<<endl;
            vector<float> row = {static_cast<float>(vir_id[i].first), static_cast<float>(vir_id[i].second), static_cast<float>(i), ans_abd[i], static_cast<float>(k)};
	    //cout<<vir_id[i].first<<" "<<vir_id[i].second<<endl;
            ANS_info.push_back(row);
            j++;
        }
    }

    return j;
}



int _Mt_Tracker::Load_id(){
     
     ifstream infile(Database.Get_Tree_Id().c_str(), ifstream::in);
     if (!infile){
                 cerr << "Error: Cannot open file : " << Database.Get_Tree_Id() << endl;
                 return 0;
                 }
     
     string buffer;
     int count = 0;
     while(getline(infile, buffer)){
                           if (buffer.size() == 0) continue;
                           Id.push_back(buffer);
                           count ++;
                           }

     infile.close();
     infile.clear();
     return count;
     }


int _Mt_Tracker::Load_order(){
    ifstream infile(Database.Get_Tree_Order().c_str(), ifstream::in);
    if (!infile){
                 cerr << "Error: Cannot open file : " << Database.Get_Tree_Order() << endl;
                 return 0;
                 }
    
     string buffer;
     int count = 0;
     while(getline(infile, buffer)){                           
                           if(buffer.size() == 0) continue;//跳过空行
                           stringstream strin(buffer);
                           int order_1 = 0;
                           int order_2 = 0;
                           int order_d = 0;
                           float dist_1 = 0;
                           float dist_2 = 0;
                           strin >> order_1 >> dist_1 >> order_2 >> dist_2 >> order_d;
                           //存储数据
                           Order_1.push_back(order_1);
                           Order_2.push_back(order_2);
                           Order_d.push_back(order_d);
                           Dist_1.push_back(dist_1);
                           Dist_2.push_back(dist_2);
                           count ++;                           
                           }

    infile.close();
    infile.clear(); 

    return count;
    }


float _Mt_Tracker::Ans_sim(float * Abd_1, vector<float> Abd_2, int j, vector<vector<float>> ANS_info){

    float Reg_1[REG_SIZE];
    float Reg_2[REG_SIZE];
    int vir[REG_SIZE];
    float total = 0;

    for(int i = 0, k = 0; i < OrderN; i++){

        int order_1 = Order_1[i];
        int order_2 = Order_2[i];
        int order_d = Order_d[i] + REG_SIZE;//Since the array subscript is non-negative, +REG_SIZE guarantees non-negative

        if(k>j-1)
            k--;
        int child_left = static_cast<int>(ANS_info[k][0]);
        int child_right = static_cast<int>(ANS_info[k][1]);

        float dist_1 = 1- Dist_1[i];//1-进化树距离
        float dist_2 = 1- Dist_2[i];
        float abd = ANS_info[k][3];



        float c1_1 = 0;//第一个菌，在第一个菌群的丰度
        float c1_2 = 0;//第一个菌，在第二个菌群的丰度

        float c2_1 = 0;//第二个菌，在第一个菌群的丰度
        float c2_2 = 0;//第二个菌，在第二个菌群的丰度

        if (order_1 >= 0){

            c1_1 = Abd_1[order_1];
            c1_2 = Abd_2[order_1];

        }
        else {
            c1_1 = Reg_1[order_1 + REG_SIZE];
            c1_2 = Reg_2[order_1 + REG_SIZE];
        }

        if (order_2 >= 0){

            c2_1 = Abd_1[order_2];
            c1_2 = Abd_2[order_2];

        }
        else {
            c2_1 = Reg_1[order_2 + REG_SIZE];
            c2_2 = Reg_2[order_2 + REG_SIZE];
        }


        float min_1 = (c1_1 < c1_2)?c1_1:c1_2;
        float min_2 = (c2_1 < c2_2)?c2_1:c2_2;

        total += min_1;
        total += min_2;

        //The abundance of common ancestors is added to the node during calculation
        if(order_1>=0&&order_2>=0){
            if(order_1==child_left&&order_2==child_right){
                Reg_2[order_d] = (c1_2 - min_1) * dist_1 + (c2_2 - min_2) * dist_2 + abd;
                k++;
            }
            else{
                Reg_2[order_d] = (c1_2 - min_1) * dist_1 + (c2_2 - min_2) * dist_2;
            }
        }
        if(order_1<0&&order_2<0){
            if(vir[order_1+ REG_SIZE]==child_left&&vir[order_2+ REG_SIZE]==child_right){
                Reg_2[order_d] = (c1_2 - min_1) * dist_1 + (c2_2 - min_2) * dist_2 + abd;
                k++;
            }
            else{
                Reg_2[order_d] = (c1_2 - min_1) * dist_1 + (c2_2 - min_2) * dist_2;
            }
        }
        if(order_1>=0&&order_2<0){
            if(order_1==child_left&&vir[order_2+ REG_SIZE]==child_right){
                Reg_2[order_d] = (c1_2 - min_1) * dist_1 + (c2_2 - min_2) * dist_2 + abd;
                k++;
            }
            else{
                Reg_2[order_d] = (c1_2 - min_1) * dist_1 + (c2_2 - min_2) * dist_2;
            }
        }
        if(order_1<0&&order_2>=0){
            if(vir[order_1+ REG_SIZE]==child_left&&order_2==child_right){
                Reg_2[order_d] = (c1_2 - min_1) * dist_1 + (c2_2 - min_2) * dist_2 + abd;
                k++;
            }
            else{
                Reg_2[order_d] = (c1_2 - min_1) * dist_1 + (c2_2 - min_2) * dist_2;
            }
        }
        Reg_1[order_d] = (c1_1 - min_1) * dist_1 + (c2_1 - min_2) * dist_2;

    }
	
    if(total!=0){
        total /= 100.0000; 
    }
    total = (total > 1.0) ? 1 : total;
    total = (total < 0.0) ? 0 : total;
    return total;
}


tuple<vector<float>, vector<float>, vector<pair<int, int>>> _Mt_Tracker::myFunction(float* Abd_1, float* Abd_2) {

    float Reg[REG_SIZE];
    float Reg_1[REG_SIZE];
    float Reg_2[REG_SIZE];
    int vir[OrderN];//Virtual code, if it is a leaf node, then it is order 1 or 2, if it is an ancestor node, then it is equal to i
    float temp_1, temp_2;

    vector<float> ANS_Leaf(OrderN + 1); // Vector array used to store leaf node abundance
    vector<float> ANS_NLeaf(OrderN + 1); // A vector storing the abundance of non-leaf nodes


    vector<pair<int, int>> ans;
    for(int i = 0; i < OrderN; i++){

        int order_1 = Order_1[i];
        int order_2 = Order_2[i];
        int order_d = Order_d[i] + REG_SIZE;
        float dist_1 = 1- Dist_1[i];
        float dist_2 = 1- Dist_2[i];

        float c1_1;
        float c1_2;

        float c2_1;
        float c2_2;

        if (order_1 >= 0){

            c1_1 = Abd_1[order_1];
            c1_2 = Abd_2[order_1];

        }
        else {
            c1_1 = Reg_1[order_1 + REG_SIZE];
            c1_2 = Reg_2[order_1 + REG_SIZE];

        }

        if (order_2 >= 0){

            c2_1 = Abd_1[order_2];
            c2_2 = Abd_2[order_2];
        }
        else {
            c2_1 = Reg_1[order_2 + REG_SIZE];
            c2_2 = Reg_2[order_2 + REG_SIZE];

        }

        

        float min_1 = (c1_1 < c1_2) ? c1_1 : c1_2;
        float min_2 = (c2_1 < c2_2)?c2_1:c2_2;


        if(order_1 >= 0){
            ANS_Leaf[order_1] = { min_1 };

        }

        if(order_2 >= 0){
            ANS_Leaf[order_2] = { min_2 };
        }


        float Nc1_1 = c1_1 - min_1;
        float Nc1_2 = c1_2 - min_1;
        float Nc2_1 = c2_1 - min_2;
        float Nc2_2 = c2_2 - min_2;
        //祖先节点丰度

        temp_1 = Nc1_1 * dist_1 + Nc2_1 * dist_2;
        temp_2 = Nc1_2 * dist_1 + Nc2_2 * dist_2;

        Reg[order_d] = (temp_1 < temp_2)?temp_1:temp_2;




        if(order_1<0 && order_2 < 0){
            ans.emplace_back(vir[order_1+ REG_SIZE],vir[order_2+ REG_SIZE]);
        }
        if(order_1>=0&&order_2<0){
            ans.emplace_back(order_1,vir[order_2+ REG_SIZE]);
        }
        if(order_1>=0&&order_2>=0){
            ans.emplace_back(order_1,order_2);
        }
        if(order_1<0&&order_2>=0){
            ans.emplace_back(vir[order_1+ REG_SIZE],order_2);
        }

        vir[order_d] = i;

        ANS_NLeaf[i] = (Reg[order_d] < 0.00001) ? 0 : Reg[order_d];

        Reg_1[order_d] = temp_1 - Reg[order_d];
        Reg_2[order_d] = temp_2 - Reg[order_d];



    }

    return make_tuple(ANS_Leaf, ANS_NLeaf, ans);
}



#endif
