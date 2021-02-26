#pragma once

#include <string>

using namespace std;

bool cmd_option_exists(char** begin, char** end, const string& option);

const char* get_cmd_option(char** begin, char** end, const string& option);
