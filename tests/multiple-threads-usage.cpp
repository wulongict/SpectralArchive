#include <iostream>
#include <thread>
#include <chrono>
#include <vector>

int main() {
  const int num_threads = std::thread::hardware_concurrency();
  
  std::cout << "Number of hardware threads: " << num_threads << std::endl;
  
  // create a vector of threads
  std::vector<std::thread> threads(num_threads);
  
  // start each thread and run it for 1 minute
  for (int i = 0; i < num_threads; ++i) {
    threads[i] = std::thread([](){
      auto start_time = std::chrono::high_resolution_clock::now();
      while (std::chrono::duration_cast<std::chrono::minutes>(std::chrono::high_resolution_clock::now() - start_time).count() < 1) {
        // do some CPU-bound work here
      }
    });
  }
  
  // wait for all threads to finish
  for (auto& thread : threads) {
    thread.join();
  }
  
  std::cout << "All threads finished." << std::endl;
  
  return 0;
}