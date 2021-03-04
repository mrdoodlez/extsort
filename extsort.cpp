/*!
 * \brief    Huge data file sorting module
 * \details  Implements sorting of binary files whos's size
 *           is bigger than internal memory amount.
 *           External sorting algorithm is used
 *
 *           https://en.wikipedia.org/wiki/External_sorting
 *
 *           Data is sorted in chunks that are swapped to disk and merged
 *           after. Source file is overwritten.
 *
 *           Sorting itself in performed using radix counting sort
 *
 *           https://www.researchgate.net/publication/283761857_Comparison_of_parallel_sorting_algorithms
 *
 *           It's proved to have best performance in case where
 *                    processors count << elements count
 *                  (no good sorting network is possible)
 *
 *           Sequential algorithm itself has complexity of k * N
 *           k - number of digits (k = 32 / 8 = 4)
 *           N - number of elemnts
 *           
 *           Given k < logN we have good performance
 *
 * \todo     Error handling
 *           || m-way merge
 *
 * \authors  Stan Raskov a.k.a mrdoodlez
 */

#include <random>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include "cmdline.h"
#include <utility>
#include <future>
#include <array>
#include <cstring>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <list>

using namespace std;

using vtype = unsigned int;

/*!
 * number of sorting passes (=k)
 */
const int pass_count = 4;

/*!
 * Number of buckets = 2 ^ (32 / 4)
 */
const int buckets_count = 256;

/*****************************************************************************/

/*!
 * \brief Class for binary partial sorted input stream
 */
class binary_istream {
  public:
    binary_istream(vtype *buffer, size_t capacity, const string& name);
    binary_istream(binary_istream&& s);
    ~binary_istream();
    void operator = (binary_istream&& s);
    vtype front() const;
    vtype evict();
    bool done() const;

  private:
    void fetch();

    vtype* buffer;
    size_t capacity;
    ifstream file;
    size_t head;
    size_t size;
    bool bdone;
    size_t fsize;
    string fname;
};

/*!
 *\ brief Class for output sorted merged stream
 */
class binary_ostream {
  public:
    binary_ostream(vtype *buffer, size_t capacity, const string& name);
    ~binary_ostream();
    void write(vtype val);
    void flush();

  private:
    vtype* buffer;
    size_t capacity;
    ofstream file;
    string fname;
    size_t head;
};

/*****************************************************************************/

/*!
 * \breief return per-thread data block size
 */
size_t get_block_size(size_t n, size_t nthreads);

/*!
 * \brief Per-thread bucket counting (prefix-sums) kernel
 */
vector<size_t> prefix_count_kernel(vtype* src, size_t n,
                                   int nthreads, int tid, int pass_num);

/*!
 * \brief Multithreaded radix sort kernel 
 */
void radix_sort_kernel(vtype* src, vtype* dst,
                       const vector<vector<size_t>>& counts, size_t n,
                       int nthreads, int tid, int pass_num);

/*!
 * \brief Radix sort caller / wrapper
 */
void radix_sort(vtype* src, vtype* back_buff, vtype** result,
                size_t n, size_t nthreads);


/*****************************************************************************/

/*!
 * \brief Random pattern data generator for self-testing
 */
vector<vtype> genrand(size_t len);

template<typename iterator>
void debug_print (iterator f, iterator l);

/*****************************************************************************/

/*!
 * \brief External sort implementation
 */
void extsort (const string& input_name, size_t memsize, size_t nthreads);

/*****************************************************************************/

/*!
 * \brief Benchmark test for different thread count
 */
void benchmark_test();

/*!
 * \brief External sort selftest
 */
void extsort_test();

/*!
 * \brief Test if input sequence is monotonically ascending
 */
template<typename iterator>
bool assert_monotonic(iterator f, iterator l);

/*****************************************************************************/

/*
 * \Print usage
 */
void usage();

/*****************************************************************************/

int main(int argc, char * argv[]) {
    if(cmd_option_exists(argv, argv + argc, "-h")) {
        usage();
        return 0;
    }
    if(cmd_option_exists(argv, argv + argc, "-b")) {
        benchmark_test();
        return 0;
    }
    
    if(cmd_option_exists(argv, argv + argc, "-e")) {
        extsort_test();
        return 0;
    }
    
    const char *fopt = get_cmd_option(argv, argv + argc, "-f");
    if (fopt == nullptr) {
        cout << "error: missing argument -v" << endl;
        usage();
        return 0;
    }
    
    const char *mopt = get_cmd_option(argv, argv + argc, "-m");
    if (mopt == nullptr) {
        cout << "error: missing argument -m" << endl;
        usage();
        return 0;
    }
    
    const char *topt = get_cmd_option(argv, argv + argc, "-t");
    if (topt == nullptr) {
        cout << "error: missing argument -t" << endl;
        usage();
        return 0;
    }

    const size_t memsize = atoi(mopt) * 1024 * 1024;
    const int nthreads = atoi(topt);

    extsort(fopt, memsize, nthreads);

    return 0;
}

/*****************************************************************************/

void extsort (const string& input_name, size_t memsize, size_t nthreads) {
    size_t data_size = memsize / sizeof(vtype) / 2;

    vtype* data      = new vtype[data_size];
    vtype* back_buff = new vtype[data_size];

    ifstream input_file(input_name, ios::binary);

    if (!input_file) {
        cout << "file " << input_name
            << " is corrupted or doesn't exist" << endl;
        return;
    }

    input_file.seekg(0, ios::end);
    size_t input_size = input_file.tellg();
    input_file.seekg(0, ios::beg);

    // sort phase

    size_t chunk_num = 0;
    while (input_size) {
        size_t chunk_size = min(data_size * sizeof(vtype), input_size);
        input_file.read((char*)data, chunk_size);
        
        auto started = chrono::high_resolution_clock::now();
        vtype* result;
        radix_sort(data, back_buff, &result,
                   chunk_size / sizeof(vtype), nthreads);
        auto done = chrono::high_resolution_clock::now();
        
        cout << "sorting chunk." << chunk_num << " of "
             << chunk_size << " bytes took "
             << chrono::duration_cast<chrono::milliseconds>(done-started).count()
             << " ms" << endl;

        ostringstream oss;
        oss << "chunk." << chunk_num;

        ofstream chunk_file(oss.str(), ios::out | ios::binary);
        chunk_file.write((char*)result, chunk_size);
        chunk_file.close();

        input_size -= chunk_size;
        chunk_num++;    
    }
    input_file.close();

    // merge phase

    size_t num_streams = chunk_num;
    size_t mem_per_stream = memsize / 2 / num_streams;
    vector<binary_istream> input_streams;  

    for (size_t stream_num = 0; stream_num < num_streams; ++stream_num) {
        ostringstream oss;
        oss << "chunk." << stream_num;    
        input_streams.push_back( {
                data + stream_num * mem_per_stream / sizeof(vtype),
                mem_per_stream / sizeof(vtype), oss.str()
        });
    }

    binary_ostream output_stream { back_buff, data_size, input_name };

    do {
        auto stream = min_element(input_streams.begin(), input_streams.end(),
        [] (const binary_istream& lhs, const binary_istream& rhs) {
            return lhs.front() < rhs.front();
        });
        output_stream.write(stream->evict());
        if (stream->done()) {
            input_streams.erase(stream);
        }

    } while(!input_streams.empty());

    output_stream.flush();

    delete[] data;
    delete[] back_buff;
}

/*****************************************************************************/

vector<vtype> genrand(size_t len) {
    static default_random_engine generator;
    uniform_int_distribution<vtype> distribution(
        numeric_limits<vtype>::min(),
        numeric_limits<vtype>::max());

    vector<vtype> v;
    while (len--) {
        v.push_back(distribution(generator));
    }
    return v;
}

template<typename iterator>
void debug_print (iterator f, iterator l) {
    for (auto i = f; i != l; ++i) {
        cout << *i << " ";
    }
    cout << endl;
}

/*****************************************************************************/

size_t get_block_size(size_t n, size_t nthreads) {
    size_t block_size = n / nthreads;
    block_size += ((n % nthreads) != 0) ? 1 : 0;
    return block_size;
}

vector<size_t> prefix_count_kernel(vtype* src, size_t n,
                                   int nthreads, int tid, int pass_num) {
    vector<size_t> count(buckets_count);
    auto block_size = get_block_size(n, nthreads);

    //LSB first bucketizer
    int offset = pass_num * 8;
    vtype mask = 0xff << offset;
    for (auto i = tid * block_size;
            (i < (tid + 1) * block_size) && (i < n); ++i) {
        int bucket = (src[i] & mask) >> offset;
        count[bucket]++;
    }
    return count;
}

void radix_sort_kernel(vtype* src, vtype* dst,
                       const vector<vector<size_t>>& counts, size_t n,
                       int nthreads, int tid, int pass_num) {
    //Find resulting prefix-sum (offset) given global counts[][]
    size_t count = 0;
    vector<size_t> prefix;
    for (int bucket = 0; bucket < buckets_count; bucket++) {
        for (int tjd = 0; tjd < nthreads; tjd++) {
            if (tjd == tid) {
                prefix.push_back(count);
            }
            count += counts[tjd][bucket];
        }
    }

    //Write element given its bucket to dst and increment prefix
    //for next element for this bucket
    auto block_size = get_block_size(n, nthreads);
    int offset = pass_num * 8;
    vtype mask = 0xff << offset;
    for (auto i = tid * block_size;
            (i < (tid + 1) * block_size) && (i < n); ++i) {
        int bucket = (src[i] & mask) >> offset;
        size_t j = prefix[bucket]++;
        dst[j] = src[i];
    }
}

void radix_sort(vtype* src, vtype* back_buff, vtype** result,
                size_t n, size_t nthreads) {
    vtype* page[2];
    int page_src = 0;
    int page_dst = 1;
    page[page_src] = src;
    page[page_dst] = back_buff;

    //4-pass stable LSB radix sort
    for (int npass = 0; npass < pass_count; npass++) {
        vector<future<vector<size_t>>> fcounts;
        //Calculate prefix-sums
        for (size_t tid = 0; tid < nthreads; tid++) {
            fcounts.push_back(async(
                                  &prefix_count_kernel, page[page_src],
                                  n, nthreads, tid, npass));
        }

        vector<vector<size_t>> counts;
        for (size_t tid = 0; tid < nthreads; tid++) {
            counts.push_back(fcounts[tid].get());
        }

        //Sort elements using prefix-sums for buckets they belong to
        vector<future<void>> fsorts;
        for (size_t tid = 0; tid < nthreads; tid++) {
            fsorts.push_back(async(
                                 &radix_sort_kernel, page[page_src],
                                 page[page_dst], counts, n,
                                 nthreads, tid, npass));
        }

        for (size_t tid = 0; tid < nthreads; tid++) {
            fsorts[tid].get();
        }

        page_src = 1 - page_src;
        page_dst = 1 - page_dst;
    }
    *result = page[page_src];
}

/*****************************************************************************/

binary_istream::binary_istream(vtype *buffer,
        size_t capacity, const string& fname) : buffer {buffer}
                                              , capacity {capacity}
                                              , file {fname, ios::binary}
                                              , head {0}
                                              , size {0}
                                              , bdone {false}
                                              , fsize {0}
                                              , fname {fname} {
    if ((buffer == 0) || (capacity == 0) || !file) {
        throw invalid_argument("binary_istream::binary_istream. invalid args");
    }

    file.seekg(0, ios::end);
    fsize = file.tellg() / sizeof(vtype);
    file.seekg(0, ios::beg);

    if (fsize == 0) {
        bdone = true;
    } else {
        fetch();
    }
}

binary_istream::binary_istream(binary_istream&& s) {
    file = move(s.file);
    fname = move(s.fname);
    buffer = s.buffer;
    capacity = s.capacity;
    head = s.head;
    size = s.size;
    fsize = s.fsize;
    bdone = s.bdone;
}

void binary_istream::operator = (binary_istream&& s) {
    file = move(s.file);
    fname = move(s.fname);
    buffer = s.buffer;
    capacity = s.capacity;
    head = s.head;
    size = s.size;
    fsize = s.fsize;
    bdone = s.bdone;

}

binary_istream::~binary_istream() {
    if (file) {
        file.close();
    }
}

vtype binary_istream::front() const {
    if (bdone) {
        return 0;
    }
    return buffer[head];    
}

vtype binary_istream::evict() {
    if (bdone) {
        return 0;
    }
    vtype ret = buffer[head++];
    if (head == size) {
        fetch();
    }
    return ret;
}

bool binary_istream::done() const {
    return bdone;
}

void binary_istream::fetch() {
    if (fsize == 0) {
        bdone = true;
        return;
    }

    size = min(fsize, capacity);
    if (!file.read((char*)buffer, size * sizeof(vtype))) {
        //TODO: handle error;
    }

    fsize -= size;
    if (fsize == 0) {
        file.close();
    }
    head = 0;
}

/*****************************************************************************/

binary_ostream::binary_ostream(vtype *buffer, size_t capacity,
        const string& fname) : buffer {buffer}
                             , capacity {capacity}
                             , file {fname, ios::binary}
                             , fname {fname}
                             , head {0} {

    if ((buffer == 0) || (capacity == 0) || !file) {
        throw invalid_argument("binary_ostream::binary_ostream. invalid args");
    }
}

binary_ostream::~binary_ostream() {
    if (file) {
        file.close();
    }
}

void binary_ostream::write(vtype val) {
    buffer[head++] = val;
    if (head == capacity) {
        flush();
    }
}

void binary_ostream::flush() {
    file.write((char*)buffer, head * sizeof(vtype));
    head = 0;
}

/*****************************************************************************/

void benchmark_test() {
    vector<size_t> test_threads = {1, 2, 4, 8, 16, 32};
    constexpr size_t data_size = 256 * 1024 * 1024;
    cout << "**** benchmark test with 256M of data ****" << endl;
    for (auto i = test_threads.begin(); i != test_threads.end(); ++i) {
        auto nthreads = *i;
        cout << "threads: " << setw(2) << nthreads << " ";
        vector<vtype> v = genrand(data_size);
        vector<vtype> back_buff(data_size);
        vtype* result;

        auto started = chrono::high_resolution_clock::now();
        radix_sort(&v[0], &back_buff[0], &result, data_size, nthreads);
        auto done = chrono::high_resolution_clock::now();

        cout << "time: "
             << chrono::duration_cast<chrono::milliseconds>(done-started).count()
             << " ms" << endl;

        //debug_print(result, result + 4);
        if (!assert_monotonic(result, result + data_size)) {
            cout << "monotonic test failed!" << endl;
        }
    }
}

void extsort_test() {
    size_t test_data_size = 1024 * 1024;
    size_t memsize = 1024 * 1024;
    int nthreads = 8;
    string fname = "test.dat";
    
    cout << "**** extsort selftest  ****" << endl;
    cout << "      data size: " << test_data_size << endl;
    cout << "   max mem size: " << memsize << " bytes" << endl;
    cout << "        threads: " << nthreads << endl;
    cout << "      data file: " << fname << endl;

    vector<vtype> v = genrand(test_data_size);
    ofstream test_file(fname, ios::binary);
    test_file.write((char*)&v[0], v.size() * sizeof(vtype));
    test_file.close();
    extsort("test.dat", memsize, nthreads);
    
    ifstream input_file(fname, ios::binary);
    if (!input_file) {
        //TODO: handle err
    }

    input_file.seekg(0, ios::end);
    size_t input_size = input_file.tellg();
    input_file.seekg(0, ios::beg);

    if (input_size / sizeof(vtype) != test_data_size) {
        cout << "sorted size assertion failed!" << endl;
        input_file.close();
        return;
    }

    input_file.read((char*)&v[0], input_size);
    input_file.close();

    if (!assert_monotonic(v.begin(), v.end())) {
        cout << "monotonic test failed!" << endl;
    }
}

template<typename iterator>
bool assert_monotonic(iterator f, iterator l) {
    auto i = f;
    auto a = *i++;
    for (; i != l; ++i) {
        if (*i < a) {
            return false;
        }
        a = *i;
    }
    return true;
}

/*****************************************************************************/

void usage() {
    cout << "usage:" << endl;
    cout << "./extsort [options]" << endl;
    cout << "options:" << endl;
    cout << "    -h                           this output" << endl;
    cout << "    -b                           perform benchmark self-test" << endl;
    cout << "    -e                           perform extsort self-test" << endl;
    cout << "    -f <input/output file name>  file to sort (in-place)" << endl;
    cout << "    -m <memory size>             available memory in MBytes" << endl;
    cout << "    -f <threads>                 thread count for radix-sort" << endl;
}
