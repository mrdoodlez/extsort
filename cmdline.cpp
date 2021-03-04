#include "cmdline.h"
#include <algorithm>

using namespace std;

bool cmd_option_exists(char** begin, char** end, const string& option) {
    return find(begin, end, option) != end;
}

const char* get_cmd_option(char** begin, char** end, const string& option) {
    char** it = find(begin, end, option);
    if (it != end && ++it != end) {
        return *it;
    }
    return nullptr;
}
