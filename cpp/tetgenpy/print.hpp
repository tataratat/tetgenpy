#pragma once

#include <iostream>
#include <stdexcept>

namespace tetgenpy {

// Adapted from:
// https://stackoverflow.com/questions/12342633/how-do-i-print-out-the-arguments-of-a-function-using-a-variadic-template

// Single Print
template<typename T>
inline void SinglePrint(T t) {
  std::cout << t << " ";
}

// Base case, no args
inline void Print() {}

// Split the parameter pack.
// We want the first argument, so we can print it.
// And the rest so we can forward it to the next call to f
template<typename T, typename... Ts>
inline void Print(T&& first, Ts&&... rest) {
  // Print first
  SinglePrint(std::forward<T>(first));
  // Forward the rest.
  Print(std::forward<Ts>(rest)...);
}

template<typename... Args>
inline void PrintInfo(Args&&... args) {
  std::cout << "TETGENPY INFO - ";
  Print(std::forward<Args>(args)...);
  std::cout << "\n";
}

/// debug printer - first argument is bool, so <on, off> is switchable.
template<typename... Args>
inline void PrintDebug(bool on, Args&&... args) {
  if (on) {
    std::cout << "TETGENPY DEBUG - ";
    Print(std::forward<Args>(args)...);
    std::cout << "\n";
  }
}

template<typename... Args>
inline void PrintWarning(Args&&... args) {
  std::cout << "TETGENPY WARNING - ";
  Print(std::forward<Args>(args)...);
  std::cout << "\n";
}

template<typename... Args>
inline void PrintAndThrowError(Args&&... args) {
  std::cout << "TETGENPY ERROR - ";
  Print(std::forward<Args>(args)...);
  std::cout << "\n";
  throw std::runtime_error("Error Occured! Abort the mission!");
}
} // namespace tetgenpy
