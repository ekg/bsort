#ifndef BSORT_HEADER_H
#define BSORT_HEADER_H

#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include <omp.h>
#include <mutex>
#include <set>
#include <dirent.h>
#include <iostream>
#include <cassert>
#include <map>
#include <vector>
#include <functional>
#include <fstream>

namespace bsort {

namespace temp_file {

/// Because the names are in a static object, we can delete them when
/// std::exit() is called.
struct Handler {
    std::set<std::string> filenames;
    std::string parent_directory;
    ~Handler() {
        // No need to lock in static destructor
        for (auto& filename : filenames) {
            std::remove(filename.c_str());
        }
        if (!parent_directory.empty()) {
            // There may be extraneous files in the directory still (like .fai files)
            auto directory = opendir(parent_directory.c_str());
            
            dirent* dp;
            while ((dp = readdir(directory)) != nullptr) {
                // For every item still in it, delete it.
                // TODO: Maybe eventually recursively delete?
                std::remove((parent_directory + "/" + dp->d_name).c_str());
            }
            closedir(directory);
            
            // Delete the directory itself
            std::remove(parent_directory.c_str());
        }
    }
};

std::string create(const std::string& base);
std::string create();
void remove(const std::string& filename);
void set_dir(const std::string& new_temp_dir);
std::string get_dir();

} // namespace temp_file

struct sort {
  int fd;
  off_t size;
  void *buffer;
};

static inline void
shellsort(unsigned char *a,
          const int n,
          const int record_size,
          const int key_size);

int cp(const char *to, const char *from);

void
radixify_parallel(const std::string& base_fn,
                  const long digit,
                  const long char_start,
                  const long char_stop,
                  const long record_size,
                  const long key_size,
                  const long stack_size,
                  const long cut_off,
                  const long switch_to_shell,
                  const long threads);

void
radixify(unsigned char *buffer,
         const long count,
         const long digit,
         const long char_start,
         const long char_stop,
         const long record_size,
         const long key_size,
         const long stack_size,
         const long cut_off,
         const long switch_to_shell);

int open_sort(char *path, struct sort *sort);

void close_sort(struct sort *sort);

}

#endif
