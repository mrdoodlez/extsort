/*!
 * \brief    Data patterns denerator
 * \authors  Stan Raskov a.k.a mrdoodlez
 */

#include <random>
#include <iostream>
#include <algorithm>
#include <fstream>
#include "cmdline.h"

using namespace std;

/*!
 * Func that prints module usage info
 */
int usage();

/*****************************************************************************/

int main(int argc, char * argv[]) {
    if(cmd_option_exists(argv, argv + argc, "-h")) {
        return usage();
    }

    const char *vopt = get_cmd_option(argv, argv + argc, "-v");
    if (vopt == nullptr) {
        cout << "error: missing argument -v" << endl;
        return usage();
    }

    const char *fopt = get_cmd_option(argv, argv + argc, "-o");
    if (fopt == nullptr) {
        cout << "error: missing argument -f" << endl;
        return usage();
    }
    
    const char *mopt = get_cmd_option(argv, argv + argc, "-m");
    if (fopt == nullptr) {
        cout << "error: missing argument -m" << endl;
        return usage();
    }

    using vtype = unsigned int;

    const size_t num_chunks = atoi(vopt);
    constexpr size_t chunk_size = 1024 * 1024 / sizeof(vtype);

    vector<vtype> chunk(chunk_size);

    default_random_engine generator;
    uniform_int_distribution<vtype> distribution(
        numeric_limits<vtype>::min(),
        numeric_limits<vtype>::max());

    vtype cnt = 0;
    vtype dcnt = chunk_size * num_chunks - 1;

    // Write output data in chunks of 1M

    if (num_chunks > 0) {
        ofstream fout(fopt, ios::out | ios::binary);
        if (*mopt == 'r') {
            for (size_t i = 0; i < num_chunks; i++) {
                generate(chunk.begin(), chunk.end(), [&distribution, &generator]() {
                    return distribution(generator);
                });
                fout.write((const char*)&chunk[0], chunk_size * sizeof(vtype));
            }
        } else if (*mopt == 'a') {
            for (size_t i = 0; i < num_chunks; i++) {
                generate(chunk.begin(), chunk.end(), [&cnt]() {
                    return cnt++;
                });
                fout.write((const char*)&chunk[0], chunk_size * sizeof(vtype));
            }
        } else if (*mopt == 'd') {
            for (size_t i = 0; i < num_chunks; i++) {
                generate(chunk.begin(), chunk.end(), [&dcnt]() {
                    return dcnt--;
                });
                fout.write((const char*)&chunk[0], chunk_size * sizeof(vtype));
            }
        } else {
            cout << "error: invlaid argument -m" << endl;
            return usage();
        }
        fout.close();
    } else {
        cout << "error: invlaid argument -v" << endl;
        return usage();
    }

    return 0;
}

/*****************************************************************************/

int usage() {
    cout << "usage:" << endl;
    cout << "./gendata [options]" << endl;
    cout << "options:" << endl;
    cout << "    -h                     this output" << endl;
    cout << "    -v <size>              output file size in MBytes" << endl;
    cout << "    -m <r|a|d>             r: random a: ascending d: descending" << endl;
    cout << "    -o <output file name>  output file name" << endl;
    return 0;
}
