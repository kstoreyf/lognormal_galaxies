#include <string>

using namespace std;

bool is_int(const string &str){
    return str.find_first_not_of("0123456789") == string::npos;
}

bool is_float(const string &str){
	signed int dec_point = str.find_first_of(".");
	if ((dec_point >= 0) and (str.find(".",dec_point+1) == string::npos)){
		   return str.find_first_not_of("0123456789.") == string::npos;
	}else{
		return 0;
	}
}

bool check_int(const string &str, const string &param_name){
	if(is_int(str)){
		return 1;
	}else{
		cout << param_name <<" should be integer!"<<endl;
		exit(1);
	}
}

bool check_float(const string &str, const string &param_name){
	if(is_float(str)){
		return 1;
	}else{
		cout << param_name <<" should be float!"<<endl;
		exit(1);
	}
}


