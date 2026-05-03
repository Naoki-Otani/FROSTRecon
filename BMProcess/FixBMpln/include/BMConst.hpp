#ifndef _INC_CONST
#define _INC_CONST

#include <array>
#include <cstddef>
#include <numeric>
#include <iostream>

// #include <spdlog/spdlog.h>

using namespace std;

/////////////////// Constants ////////////////////

// number of MCR
static const std::size_t NUM_MCR = 8;
// number of slot in each MCR
static const std::size_t NUM_MCR_SLOT = 8;
// number of total slot
static const std::size_t NUM_TOTAL_SLOT = NUM_MCR * NUM_MCR_SLOT;
// number of FEB
static const std::size_t NUM_FEB = 44;
// number of chip per FEB
static const std::size_t NUM_CHIP_PER_FEB = 3;
// number of channel per chip
static const std::size_t NUM_CHANNEL_PER_CHIP = 32;
// number of FEB in each MCR
static constexpr std::array<std::size_t, NUM_MCR> NUM_FEB_PER_MCR = {5, 6, 6, 6, 5, 5, 6, 5};
// FEB inde -> FEB ID
static constexpr std::array<int, NUM_FEB> FEB_ID = {
   0,  1,  2,  3,  4,        // MCR 0
   8,  9, 10, 11, 12, 13,    // MCR 1
  16, 17, 18, 19, 20, 21,    // MCR 2
  24, 25, 26, 27, 28, 29,    // MCR 3
  32, 33, 34, 35, 36,        // MCR 4
  40, 41, 42, 43, 44,        // MCR 5
  48, 49, 50, 51, 52, 53,    // MCR 6
  56, 57, 58, 59, 60         // MCR 7
};
// number for invalid feb index
constexpr std::size_t INVALID_FEB_INDEX = 999;
// FEB ID -> number of channels
static constexpr std::array<std::size_t, NUM_MCR * NUM_MCR_SLOT> FEB_NUM_CHANNEL = {
  95, 95, 95, 95, 95,  0,  0,  0,   //  0– 7
  95, 95, 95, 95, 95, 96,  0,  0,   //  8–15
  95, 95, 95, 95, 96, 14,  0,  0,   // 16–23
  95, 95, 95, 95, 96, 14,  0,  0,   // 24–31
  95, 95, 95, 95, 96,  0,  0,  0,   // 32–39
  95, 95, 95, 95, 96,  0,  0,  0,   // 40–47
  95, 95, 95, 95, 96, 95,  0,  0,   // 48–55
  95, 95, 95, 95, 95,  0,  0,  0    // 56–63
};
// FEB ID for YASU tracker
static constexpr std::array<int, 2> YASU_FEB_ID = {21, 29};

// number of horizontal scintillator bar in one plane
static const std::size_t NUM_HOR_BAR = 95;
// number of vertical scintillator bar in one plane
static const std::size_t NUM_VER_BAR = 16;
// number of scintillator bar in one YASU plane
static const std::size_t NUM_YASU_BAR = 7;

// number of planes of BabyMIND
static const std::size_t NUM_BM_PLANE = 18;
// number of planes of YASU tracker
static const std::size_t NUM_YASU_PLANE = 2;


/////////////////// Converter ////////////////////

// MCR number and Slot number → FEB ID
constexpr int MakeFebID(std::size_t mcr, std::size_t slot) {
  return mcr * NUM_MCR_SLOT + slot;
}
// Check if FEB with given FEB ID exist
constexpr bool IsValidFebID(int feb_id){
  for (int id : FEB_ID)
    if (id == feb_id) return true;

  return false;
}
// Converter from FEB ID to index
constexpr std::size_t FebIDToIndex(int feb_id) {
  for (std::size_t index = 0; index < FEB_ID.size(); ++index)
    if (FEB_ID[index] == feb_id) return index;

  // spdlog::warn("Invalid FEB ID : {}", feb_id);
  // std::cout << "Invalid FEB ID : " << feb_id << std::endl;
  return INVALID_FEB_INDEX;
}

#endif
