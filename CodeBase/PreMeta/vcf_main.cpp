// Date: March 2015
// Author: paulbunn@email.unc.edu (Paul Bunn)

#include <ctime>     // For printing out current time
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

bool ReadFileViaGetline(const string& filename,
                        const bool do_processing,
                        int64_t* num_lines) {
  string line;
  time_t current_time = time(nullptr);
  cout << asctime(localtime(&current_time))
       << "ReadFileViaGetline. Opening file..." << endl;
  ifstream myfile (filename);
  current_time = time(nullptr);
  cout << asctime(localtime(&current_time))
       << "ReadFileViaGetline. File open. Reading..." << endl;
  if (myfile.is_open()) {
    while (getline(myfile, line)) {
      if (do_processing) {
        (*num_lines)++;
        //cout << line << '\n';
      }
    }
    current_time = time(nullptr);
    cout << asctime(localtime(&current_time))
         << "ReadFileViaGetline. Done Reading " << *num_lines
         << " lines. Closing file..." << endl;
    myfile.close();
    return true;
  }

  else {
    current_time = time(nullptr);
    cout << asctime(localtime(&current_time))
         << "ReadFileViaGetline. Unable to open file." << endl;
    return false;
  }
}

bool ReadFileViaVectorBuffer(const string& filename,
                             const bool do_processing_iteration,
                             int64_t* num_lines) {
  //const int64_t max_buffer_size = 214748364; // .2 GiB
  const int64_t max_buffer_size = 2147483648; // 2 GiB

  time_t current_time = time(nullptr);
  cout << asctime(localtime(&current_time))
       << "ReadFileViaVectorBuffer. Opening file..." << endl;
  // Open file at end (ios::ate) so we can first grab its size.
  ifstream file (filename, ios::in|ios::binary|ios::ate);
  if (file.is_open()) {
    int64_t file_size = file.tellg();
    const int64_t block_size = min(max_buffer_size, file_size);
    cout << "Total file size: " << file_size << ", block_size: "
         << block_size << endl;

    // Return to beginning of file.
    file.seekg(0, ios::beg);
    int64_t iteration = 1;
    int64_t last_endline_character_pos = 0;
    int64_t prev_start_pos = 0;
    while (true) {
      // Print time.
      current_time = time(nullptr);
      cout << asctime(localtime(&current_time))
           << "ReadFileViaVectorBuffer. Reading Block " << iteration
           << " to memory. Previous start pos: " << prev_start_pos
           << ", Current start pos: " << last_endline_character_pos
           << "..." << endl;

      // Set start position based on last '\n' character read from previous block.
      if (file.tellg() != 0) {
        file.seekg(last_endline_character_pos);
      }
      prev_start_pos = file.tellg();

      const int64_t remaining_bytes = file_size - file.tellg();
      // Check if we're done reading file.
      if (remaining_bytes <= 0) {
        cout << "Aborting while loop: non-positive remaining_bytes: "
             << remaining_bytes << endl;
        break;
      }

      // Reserve space for memblock: either block_size OR remaining_bytes,
      // where we use block_size unless this is the final
      // block and there are fewer than block_size bytes left.
      int64_t current_block_size = min(remaining_bytes, block_size);
      vector<char> memblock;
      memblock.reserve(current_block_size);

      // Read block.
      memblock.assign(istreambuf_iterator<char>(file), istreambuf_iterator<char>());

      // Process block.
      current_time = time(nullptr);
      cout << asctime(localtime(&current_time))
           << "ReadFileViaVectorBuffer. Done Reading Block " << iteration
           << "Now Processing it..." << endl;
      if (do_processing_iteration) {
        for (int64_t i = 0; i < memblock.size(); i++) {
          if (memblock[i] == '\n') {
            last_endline_character_pos = file.tellg();
            (*num_lines)++;
          } else if (memblock[i] == '\0') {
            cout << "Aborting processing early (" << i << " vs. "
                 << memblock.size() << ")" << endl;
            break;
          }
        }
      }
      current_time = time(nullptr);
      cout << asctime(localtime(&current_time))
           << "ReadFileViaVectorBuffer. Done Processing Block " << iteration
           << ". Proceeding to next Block..." << endl;
      iteration++;
    }

    current_time = time(nullptr);
    cout << asctime(localtime(&current_time))
         << "ReadFileViaVectorBuffer. Done Processing " << *num_lines
         << " lines. Closing file..." << endl;
    file.close();
    current_time = time(nullptr);
    cout << asctime(localtime(&current_time))
         << "ReadFileViaVectorBuffer. File closed. Deleting buffer..." << endl;
    return true;
  }
  else {
    current_time = time(nullptr);
    cout << asctime(localtime(&current_time))
         << "ReadFileViaVectorBuffer. Unable to open file." << endl;
    return false;
  }
}

// TODO(PHB): Make this do something useful.
void ProcessLine(const string& line) {
  // Remove trailing CR (Carriage Return, ascii value = 13) character,
  // if necessary.
  if (line[line.length() - 1] == 13) {
    //cout << "|" << line.substr(0, line.length() - 1) << "|" << endl;
  } else {
    //cout << "Line|" << line << "|" << endl;
  }
}

bool ReadFileViaCharBuffer(const string& filename,
                           const bool do_processing_stream_getline,
                           const bool do_processing_stream_manual,
                           int64_t* num_lines) {
  //const int64_t max_buffer_size = 614748364; // .6 GiB
  const int64_t max_buffer_size = 2147483648; // 2 GiB

  time_t current_time = time(nullptr);
  cout << asctime(localtime(&current_time))
       << "ReadFileViaCharBuffer. Opening file..." << endl;
  // Open file at end (ios::ate) so we can first grab its size.
  ifstream file (filename, ios::in|ios::binary|ios::ate);
  if (file.is_open()) {
    int64_t file_size = file.tellg();
    const int64_t block_size = min(max_buffer_size, file_size);
    char* memblock = new char[block_size];
    cout << "Total file size: " << file_size << ", block_size: "
         << block_size << endl;

    // Return to beginning of file.
    file.seekg(0, ios::beg);
    int64_t iteration = 1;
    int64_t next_block_start_pos = file.tellg();
    int64_t prev_start_pos = -1;
    string final_line;
    while (true) {
      // Print time.
      current_time = time(nullptr);
      cout << asctime(localtime(&current_time))
           << "ReadFileViaCharBuffer. Reading Block " << iteration
           << " to memory. Previous start pos: " << prev_start_pos
           << ", Current start pos: " << next_block_start_pos
           << "..." << endl;

      // Sanity check we've made progress since the last block was read.
      int64_t current_file_pos = file.tellg();
      if (prev_start_pos == next_block_start_pos) {
        if (prev_start_pos >= 0) {
          // We enter here when either:
          //  a) The last line from the previous block was so long that there
          //     wasn't a single line-break in the entire block; OR
          //  b) We reached the last line of the file.
          // In either case, we store the last line, and abort.
          if (!final_line.empty()) {
            (*num_lines)++;
            ProcessLine(final_line);
          }
          if (current_file_pos + final_line.length() == file_size) {
            // Case (b): This is the last line of the file.
            return true;
          } else {
            // Case (a): This line is too long (fills up the entire buffer before
            // reaching the NewLine character.
            cout << "ERROR: Line starting at position " << current_file_pos
                 << " is too long. Aborting." << endl;
            return false;
          }
        }
      } else if (current_file_pos != next_block_start_pos) {
        file.seekg(next_block_start_pos);
      }

      prev_start_pos = file.tellg();
      const int64_t remaining_bytes = file_size - prev_start_pos;
      // Check if we're done reading file.
      if (remaining_bytes <= 0) {
        cout << "Aborting while loop: non-positive remaining_bytes: "
             << remaining_bytes << endl;
        break;
      }

      // Resize memblock, if this is the final block and there are fewer
      // than block_size bytes left.
      int64_t current_block_size = min(remaining_bytes, block_size);
      if (current_block_size < block_size) {
        delete[] memblock;
        memblock = new char[current_block_size];
      }

      // Read block.
      file.read(memblock, current_block_size);

      // Process block.
      current_time = time(nullptr);
      cout << asctime(localtime(&current_time))
           << "ReadFileViaCharBuffer. Done Reading Block " << iteration
           << ". Now Processing it..." << endl;
      if (do_processing_stream_getline) {
        stringstream ss;
        ss.rdbuf()->pubsetbuf(memblock, current_block_size);

        // Determine if this block ends right at the end of line, as
        // determined by the New Line character '10'.
        ss.seekg(-1, ios::end);
        int end_char = ss.peek();
        bool last_char_is_newline = end_char == 10;
        ss.seekg(0, ios::beg);

        // Process Lines.
        string str;
        int64_t current_stream_pos;
        bool at_least_one_line = false;
        bool is_first_line_of_block = true;
        while(getline(ss, str)) {
          if (is_first_line_of_block) {
            is_first_line_of_block = false;
            str = final_line + str;
          }
          current_stream_pos = ss.tellg();
          if (current_stream_pos == -1) {
            final_line = str;
          } else {
            // Not the last line of the block. Update num_lines.
            (*num_lines)++;
            at_least_one_line = true;
            ProcessLine(str);
          }
        }

        // If last character in the block was NL character, then final_line
        // will not get reset from its value from the previous block. Reset
        // it now.
        if (last_char_is_newline) {
          final_line = "";
        }

        // Set next block start position.
        if (at_least_one_line) {
          next_block_start_pos = file.tellg();
        }
      } else if (do_processing_stream_manual) {
        stringstream ss;
        ss.rdbuf()->pubsetbuf(memblock, current_block_size);
        // Determine if this block ends right on an end of line: these
        // have form \13 \10 (13 = Carriage Return, 10 = New Line).
        // If not, then we'll drop the last line we read from this
        // block, and have the next block re-fetch that line. If so,
        // we'll keep the last line. Some extra code is required to
        // handle the case that we're exactly between the \13 and the
        // \10 (i.e. the block ends with \13 character).
        int num_carriage_returns = 0;
        int char_int = ss.get();
        while (char_int >= 0) {
          if (char_int == 13) {
            num_carriage_returns++;
          } else if (char_int == 10) {
            (*num_lines)++;
          }
          char_int = ss.get();
        }

        next_block_start_pos = file.tellg();

        // Very last line of file wasn't counted, since it didn't end
        // in NewLine character. Count it now.
        if (next_block_start_pos == -1) {
          (*num_lines)++;
        }
      }

      // Record total processing time for this block.
      current_time = time(nullptr);
      cout << asctime(localtime(&current_time))
           << "ReadFileViaCharBuffer. Done Processing Block " << iteration
           << ". Proceeding to next Block..." << endl;
      iteration++;
    }

    // Record Final time.
    current_time = time(nullptr);
    cout << asctime(localtime(&current_time))
         << "ReadFileViaCharBuffer. Done Processing " << *num_lines
         << " lines. Closing file..." << endl;
    file.close();
    current_time = time(nullptr);
    cout << asctime(localtime(&current_time))
         << "ReadFileViaCharBuffer. File closed. Deleting buffer..." << endl;
    delete[] memblock;
    return true;
  }
  else {
    current_time = time(nullptr);
    cout << asctime(localtime(&current_time))
         << "ReadFileViaCharBuffer. Unable to open file." << endl;
    return false;
  }
}

int main(int argc, char* argv[]) {
  // The 244G .vcf file
  const string& filename = "G:/projects/lin/ESP6800/vcf/combined.vcf";
  // A 2G sample:
  //const string& filename = "G:/projects/lin/2012-05-10/6900/vcf_pheno.txt";
  // A 24 byte sample:
  //const string& filename = "H:/Documents/foo.trash";
  int64_t num_lines = 0;
  bool do_getline = false;
  bool do_char_buffer_stream_getline = true;
  bool do_char_buffer_stream_manual = false;
  bool do_vector_buffer_iteration = false;
  time_t current_time;

  // Read file using GetLine.
  if (do_getline) {
    current_time = time(nullptr);
    cout << endl << asctime(localtime(&current_time)) << "Reading via getline..."
         << endl << endl;
    if (!ReadFileViaGetline(filename, true, &num_lines)) {
      cout << "Failed ReadFileViaGetline." << endl;
    }
    current_time = time(nullptr);
    cout << endl << asctime(localtime(&current_time))
         << "Done reading via getline." << endl;
  }

  // Read file using char buffer via stream.
  if (do_char_buffer_stream_getline) {
    num_lines = 0;
    cout << endl << "Reading File via char buffer via stream..." << endl;
    if (!ReadFileViaCharBuffer(filename, true, false, &num_lines)) {
      cout << "Failed ReadFileViaCharBuffer via stream." << endl;
    }
    current_time = time(nullptr);
    cout << endl << asctime(localtime(&current_time))
         << "Done reading char buffer via stream." << endl;
  }

  // Read file using char buffer via stream.
  if (do_char_buffer_stream_manual) {
    num_lines = 0;
    cout << endl << "Reading File via char buffer via stream..." << endl;
    if (!ReadFileViaCharBuffer(filename, false, true, &num_lines)) {
      cout << "Failed ReadFileViaCharBuffer via stream." << endl;
    }
    current_time = time(nullptr);
    cout << endl << asctime(localtime(&current_time))
         << "Done reading char buffer via stream." << endl;
  }

  // Read file using vector buffer via iterating.
  if (do_vector_buffer_iteration) {
    num_lines = 0;
    cout << endl << "Reading File via vector buffer via iteration..." << endl;
    if (!ReadFileViaVectorBuffer(filename, true, &num_lines)) {
      cout << "Failed ReadFileViaVectorBuffer via iteration." << endl;
    }
    current_time = time(nullptr);
    cout << endl << asctime(localtime(&current_time))
         << "Done reading via vector buffer via iteration." << endl;
  }

  return 0;
}
