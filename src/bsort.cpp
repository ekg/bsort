#include "bsort.hpp"

namespace bsort {


namespace temp_file {

// We use this to make the API thread-safe
std::recursive_mutex monitor;
std::string temp_dir;
Handler handler;

std::string create(const std::string& base) {
    std::lock_guard<std::recursive_mutex> lock(monitor);

    if (handler.parent_directory.empty()) {
        // Make a parent directory for our temp files
        std::string tmpdirname = get_dir() + "/bsort-XXXXXX";
        auto got = mkdtemp(&tmpdirname[0]);
        if (got != nullptr) {
            // Save the directory we got
            handler.parent_directory = got;
        } else {
            std::cerr << "[bsort]: couldn't create temp directory: " << tmpdirname << std::endl;
            exit(1);
        }
    }

    std::string tmpname = handler.parent_directory + "/" + base + "XXXXXX";
    // hack to use mkstemp to get us a safe temporary file name
    int fd = mkstemp(&tmpname[0]);
    if(fd != -1) {
        // we don't leave it open; we are assumed to open it again externally
        close(fd);
    } else {
        std::cerr << "[bsort]: couldn't create temp file on base "
                  << base << " : " << tmpname << std::endl;
        exit(1);
    }
    handler.filenames.insert(tmpname);
    return tmpname;
}

std::string create() {
    // No need to lock as we call this thing that locks
    return create("bsort-");
}

void remove(const std::string& filename) {
    std::lock_guard<std::recursive_mutex> lock(monitor);
    
    std::remove(filename.c_str());
    handler.filenames.erase(filename);
}

void set_dir(const std::string& new_temp_dir) {
    std::lock_guard<std::recursive_mutex> lock(monitor);
    
    temp_dir = new_temp_dir;
}

std::string get_dir() {
    std::lock_guard<std::recursive_mutex> lock(monitor);

    // Get the default temp dir from environment variables.
    if (temp_dir.empty()) {
        const char* system_temp_dir = nullptr;
        for(const char* var_name : {"TMPDIR", "TMP", "TEMP", "TEMPDIR", "USERPROFILE"}) {
            if (system_temp_dir == nullptr) {
                system_temp_dir = getenv(var_name);
            }
        }
        temp_dir = (system_temp_dir == nullptr ? "/tmp" : system_temp_dir);
    }

    return temp_dir;
}

} // namespace temp_file

static inline void
shellsort(unsigned char *a,
          const int n,
          const int record_size,
          const int key_size) {
  int i, j;
  char temp[record_size];

  for (i=3; i < n; i++) {
    memcpy(&temp, &a[i * record_size], record_size);
    for(j=i; j>=3 && memcmp(a+(j-3)*record_size, &temp, key_size) >0; j -= 3) {
      memcpy(a+j*record_size, a+(j-3)*record_size, record_size);
    }
    memcpy(a+j*record_size, &temp, record_size);
  }

  for (i=1; i < n; i++) {
    memcpy(&temp, &a[i*record_size], record_size);
    for(j=i; j>=1 && memcmp(a+(j-1)*record_size, &temp, key_size) >0; j -= 1) {
      memcpy(a+j*record_size, a+(j-1)*record_size, record_size);
    }
    memcpy(a+j*record_size, &temp, record_size);
  }
}

int cp(const std::string& from, const std::string& to) {
    std::ofstream out(to, std::ios_base::binary);
    std::ifstream in(from, std::ios_base::binary);
    out << in.rdbuf();
    in.close();
    out.close();
}

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
                  const long threads) {

    temp_file::temp_dir = "./";
    // copy the base file to a new tempfile
    std::string temp_fn = temp_file::create(std::string(base_fn));
    cp(base_fn, temp_fn);
    //std::cerr << "temp file " << temp_fn << std::endl;
    // open the tempfile
    sort temp_s;
    if (-1==bsort::open_sort((char*)temp_fn.c_str(), &temp_s)) {
        assert(false);
    }
  
    // compute starting positions of buffers
    size_t chunk_max = temp_s.size/threads + temp_s.size%threads;
    //std::cerr << "Chunk max " << chunk_max << std::endl;
    unsigned char* buffer[threads];
    size_t chunk_size[threads];
    size_t start = 0;
    for (long i = 0; i < threads; ++i) {
        buffer[i] = (unsigned char*)(temp_s.buffer) + start;
        chunk_size[i] = (start + chunk_max > temp_s.size ? temp_s.size - start : chunk_max);
        start += chunk_size[i];
    }
    /*
    for (long i = 0; i < threads; ++i) {
        std::cerr << " chunk start and size " << &buffer[i] << " " << chunk_size[i] << std::endl;
    }
    */
    //std::cerr << "start " << start << " temp_s.size " << temp_s.size << std::endl;
    assert(start == temp_s.size);

    //omp_set_nested(true);
    
#pragma omp parallel for schedule(dynamic)
    for (long i = 0; i < threads; ++i) {
#pragma omp task
        radixify(buffer[i], //(unsigned char*)temp_s.buffer+offset[i], // todo get buffer pointer for this chunk
                 chunk_size[i]/record_size, // todo get count for this chunk
                 digit,
                 char_start,
                 char_stop,
                 record_size,
                 key_size,
                 stack_size,
                 cut_off,
                 switch_to_shell);

    }

    // sorts for each chunk are done

    // open the base file
    sort base_s;
    if (-1==bsort::open_sort((char*)base_fn.c_str(), &base_s)) {
        assert(false);
    }

    // map from key to offset/reader pointers
    struct cmp_char_ptr {
        bool operator()(const std::pair<unsigned char*, size_t>& a, const std::pair<unsigned char*, size_t>& b) const {
            return memcmp(a.first, b.first, a.second) < 0;
        }
    };
    std::map<std::pair<unsigned char*, size_t>, std::vector<std::pair<size_t, size_t> >, cmp_char_ptr> priority_q;

    // fill the priority queue
    for (long i = 0; i < threads; ++i) {
        priority_q[std::make_pair(&buffer[i][0], key_size)].push_back(std::make_pair(i, 0));
    }
    size_t write_at = 0;
    while (!priority_q.empty()) {
        auto b = priority_q.begin();
        std::vector<std::pair<size_t, size_t> > ps = b->second;
        for (auto& p : ps) {
            memcpy((unsigned char*)base_s.buffer+write_at,
                   (unsigned char*)buffer[p.first]+p.second, record_size);
            write_at += record_size;
        }
        priority_q.erase(priority_q.begin());
        for (auto& p : ps) {
            auto new_off = p.second + record_size;
            if (new_off + record_size <= chunk_size[p.first]) {
                unsigned char* c = (unsigned char*)buffer[p.first]+new_off;
                priority_q[std::make_pair(c, key_size)].push_back(std::make_pair(p.first, new_off));
            }
        }
    }

    // remove our temp file
    std::remove(temp_fn.c_str());
    close_sort(&base_s);

}

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
         const long switch_to_shell) {
  long counts[char_stop+1];
  long offsets[char_stop+1];
  long starts[char_stop+1];
  long ends[char_stop+1];
  long offset=0;
  unsigned char temp[record_size];
  long target, x, a, b;
  long stack[stack_size];
  long stack_pointer;
  long last_position, last_value, next_value;

  for (x=char_start; x<=char_stop; x++) {
    counts[x] = 0;
    offsets[x] = 0;
  }

  // Compute starting positions
  for (x=0; x<count; x++) {
    long c = buffer[x*record_size + digit];
    counts[c] += 1;
  }

  // Compute offsets
  offset = 0;
  for(x=char_start; x<=char_stop; x++) {
    offsets[x] = offset;
    starts[x] = offsets[x];
    offset += counts[x];
  }

  for(x=char_start; x<char_stop; x++) {
    ends[x] = offsets[x+1];
  }
  ends[char_stop] = count;

  for(x=char_start; x<=char_stop; x++) {
    while (offsets[x] < ends[x]) {

      if (buffer[offsets[x] * record_size + digit] == x) {
        offsets[x] += 1;
      } else {
        stack_pointer=0;
        stack[stack_pointer] = offsets[x];
        stack_pointer += 1;
        target = buffer[offsets[x] * record_size + digit];
        while( target != x && stack_pointer < stack_size ) {
          stack[stack_pointer] = offsets[target];
          offsets[target] += 1;
          target = buffer[stack[stack_pointer] * record_size + digit];
          stack_pointer++;
        };
        if (stack_pointer != stack_size) {
          offsets[x] += 1;
        }
        stack_pointer--;
        memcpy(&temp, &buffer[stack[stack_pointer] * record_size], record_size);
        while (stack_pointer) {
          memcpy(&buffer[stack[stack_pointer] * record_size], &buffer[stack[stack_pointer-1] * record_size], record_size);
          stack_pointer--;
        }
        memcpy(&buffer[stack[0] * record_size], &temp, record_size);
      }
    }
  }

#pragma omp parallel
#pragma omp single nowait
  {
    if (digit < cut_off) {
      for(x=char_start; x<=char_stop; x++) {
        if ( ends[x] - starts[x] > switch_to_shell) {
#pragma omp task shared(buffer)
          radixify(&buffer[starts[x] * record_size],
                   ends[x] - starts[x],
                   digit+1,
                   char_start,
                   char_stop,
                   record_size,
                   key_size,
                   stack_size,
                   cut_off,
                   switch_to_shell);
        } else {
          if (ends[x] - starts[x] <= 1) continue;
#pragma omp task shared(buffer)
          shellsort(&buffer[starts[x] * record_size], ends[x] - starts[x], record_size, key_size);
        }
      }
    } else {
      for(x=char_start; x<=char_stop; x++) {
        if (ends[x] - starts[x] > 1) {
#pragma omp task shared(buffer)
          shellsort(&buffer[starts[x] * record_size], ends[x] - starts[x], record_size, key_size);
        }
      }
    }
  }
}


int open_sort(char *path, struct sort *sort) {
  void *buffer = NULL;

  int fd = open(path, O_RDWR);
  if (fd == -1)
    goto error;

  struct stat stats;
  if (-1 == fstat(fd, &stats))
    goto error;
  if (!(buffer = mmap(NULL,
                      stats.st_size,
                      PROT_READ | PROT_WRITE,
                      MAP_SHARED,
                      fd,
                      0
                     )))
    goto error;


  madvise(buffer, stats.st_size, POSIX_MADV_WILLNEED | POSIX_MADV_SEQUENTIAL);

  sort->buffer = buffer;
  sort->size = stats.st_size;
  sort->fd = fd;
  return 0;

error:
  perror(path);
  if (buffer)
    munmap(buffer, stats.st_size);
  if (fd != -1)
    close(fd);
  sort->buffer = 0;
  sort->fd = 0;
  return -1;
}


void close_sort(struct sort *sort) {
  if (sort->buffer) {
    munmap(sort->buffer, sort->size);
    sort->buffer = 0;
    sort->size = 0;
  }

  if (sort->fd) {
    close(sort->fd);
    sort->fd = 0;
  }
}

}
