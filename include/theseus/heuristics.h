/*
 *                             The MIT License
 *
 * Copyright (c) 2024 by Albert Jimenez-Blanco
 *
 * This file is part of #################### Theseus Library ####################.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 */

#pragma once

#include <iostream>
#include <vector>

#include "theseus/alignment.h"

/**
 * Class containing the data related to heuristics
 *
 */

namespace theseus {

    class Heuristics {
    public:
        /**
         * @brief Construct a new Heuristics object with parameters
         *
         */
        Heuristics(
            // int  s_min,
            int  lag_threshold,
            int  lag_density_threshold,
            double density_threshold,
            int  max_steps
        ) {
            // Used heuristics
            _lag_behind          = (lag_threshold != -1);
            _lag_density         = (lag_density_threshold != -1);
            _advancement_density = (density_threshold != -1);
            _steps_limit         = max_steps != -1;

            // Heuristics' parameters
            _lag_threshold       = lag_threshold;
            _density_threshold   = density_threshold;
            _max_steps           = max_steps;
            _min_offset_increase = lag_density_threshold;

            // Max offset default
            _max_offset          = 0;
        }


        // INITIALIZER //
        void new_alignment(int scope_size, int mismatch, int gape) {
            // General
            _max_offset = 0;

            // Initialize the vector of the last max offsets to 0
            _K = 10*scope_size;
            _scores_to_affect = 0;
            _advanced_enough = false;
            if (_last_max_offsets.size()==0) {
                int num_mismatches = (10*scope_size)/mismatch;
                int num_gaps = (10*scope_size)/(gape);
                _min_offset_increase += + std::max(num_mismatches, num_gaps);
                _last_max_offsets.resize(_K, 0);
            }
            else {
                std::fill(_last_max_offsets.begin(), _last_max_offsets.end(), 0);
            }
        }


        // ------------------------ HEURISTICS ---------------------------------

        // LOCAL heuristics: Heuristics that are checked on all wavefront cells
        // TODO: How do the two local heuristics coexist?
        bool check_local_heuristics(int &curr_offset) {
            // Check lag density
            if (_lag_density && check_lag_density(curr_offset)) {
                return true;
            }
            else if (_lag_behind && check_lag_behind(curr_offset)) {
                return true;
            }
            else return false;
        }

        /**
         * @brief Lag behind heuristic
         *
         * Let "max_offset" indicate the maximum offset reached so far by Theseus
         * and let "curr_offset" indicate the offset of a currently analysed Theseus'
         * cell. Also let "lag_threshold" indicate a maximum allowed offset difference
         * between "max_offset" and "curr_offset".
         * -----------
         * Definition: The lag behing heuristic checks, for a given cell with offset
         * curr_offset whether
         *                      max_offset - curr_offset > lag_threshold
         * That is, whether the current cell if too far behind from the currently
         * most promising cell (the one with offset max_offset).
         * -----------
         * We can use this heuristic to prune unpromising cells and paths.
         */
        inline bool check_lag_behind(int &curr_offset) {
            if (!_lag_behind) return false;
            _max_offset = (curr_offset > _max_offset)? curr_offset : _max_offset;
            return (_max_offset - curr_offset > _lag_threshold);
        }

        /**
         * @brief  Density lag heuristic.
         *
         * This heuristic combines is a refined version of the lag behind heuristic,
         * that only prunes cells when they are too far behind AND the most promising
         * path has advanced enough in the last K steps. That is, you want to avoid
         * pruning cells when the most promising path has not advanced in the last K
         * steps, as it may be stuck in a zone with a high error rate in the alignment.
         * -----------
         * Definition: We define the density lag heuristic as a heuristic that checks,
         * for a given cell with offset curr_offset if two conditions hold:
         *      1. max_offset - curr_offset > lag_threshold
         *      2. The most promising path has advanced enough in the last K steps
         *         (that is, the increase in the maximum offset in the last K steps
         *         is greater than a given threshold, indicating a "promising advancement
         *         density").
         *            max_offset - _last_max_offsets[score%K] > min_increase
         * -----------
         * We can use this heuristic to prune unpromising cells and paths, while avoiding
         * pruning cells when there is no clear primising path.
         */
        inline bool check_lag_density(int &curr_offset) {
            if (!_lag_density) return false;
            // Update the max value
            _max_offset = (curr_offset > _max_offset)? curr_offset : _max_offset;
            // Condition 1
            bool condition_1 = _max_offset - curr_offset > _min_offset_increase;
            // Condition 2 is _advanced_enough value (updated once per score)
            return (condition_1 && _advanced_enough);
        }


        // GLOBAL heuristics: Heuristics checked once per score
        int check_global_heuristics(int score) {
            // Update density lag heuristic
            if (_lag_density) update_density(score);
            // Check advancement density heuristic
            if (_advancement_density) return check_advancement_density(score);
            // Check step limit heuristic
            if (_steps_limit) return is_over_step_limit(score);
            return THESEUS_STATUS_OK;
        }

        /**
         * @brief Advancement density
         *
         * Let "max_offset" indicate the maximum offset reached so far by Theseus
         * and let "density_threshold" indicate the minimum advancement density to
         * continue the alignment.
         * -----------
         * Definition: We define the advancement density heuristic as a heuristic
         * that checks whether the current offset/score ratio is lower than a given
         * threshold. That is, we check
         *                      max_offset/score < density_threshold
         * We further complement this by checking if the score is greater than a
         * minimum score "s_min", to avoid stopping alignment with little context.
         * -----------
         * This heuristic serves as a global stopping criterion: if the heuristic
         * fails, Theseus stops alignment and returns an incomplete alignment.
         */
        inline int check_advancement_density(int &score) {
            // Check if the heuristic is active
            if (!_advancement_density) return THESEUS_STATUS_OK;
            // Check minimum score
            if (score <= _s_min) return THESEUS_STATUS_OK;
            // Return heuristic fail or success
            return (static_cast<double>(_max_offset)/score < _density_threshold) ? THESEUS_STATUS_END_UNREACHABLE : THESEUS_STATUS_OK;
        }

        /**
         * @brief Update density for the lag density heuristic
         *
         * @param score
         */
        void update_density(int &score) {
            // Check if density is high  enough
            if (!_lag_density) return;
            // Before updating, _last_max_offsets[score%_K] contains the max
            // offset reached K steps ago.
            int offset_inc = _max_offset - _last_max_offsets[score%_K];
            if (offset_inc > _min_offset_increase) {
                _advanced_enough = true;
                // _scores_to_affect = _K/10; // Once we have a promising advancement, the heuristic affects the next K scores
            }
            else {
                _advanced_enough = false;
            }
            // _scores_to_affect -= 1;
            // if (_scores_to_affect <= 0) {
            //     _advanced_enough = false; // If we have affected K scores, we stop affecting scores until we find another promising advancement
            // }
            // Update the vcalue for the max offset in the current score
            _last_max_offsets[score%_K] = _max_offset;
            // std::cout << "advanced_enough: " << _advanced_enough << " offset_inc: " << offset_inc << " max_offset: " << _max_offset << std::endl;
        }

        /**
         * @brief Over steps limit heuristic
         *
         * Check if we are out of the step limit (the two sequences are too
         * similar)
         *
         */
        int is_over_step_limit(int &score) {
            if (!_steps_limit) return THESEUS_STATUS_OK;
            else return score > _max_steps ? THESEUS_STATUS_MAX_STEPS_REACHED : THESEUS_STATUS_OK;
        }


        //---------------------------- ACCESSORS -------------------------------

        /**
         * @brief Return whether the lag behind heuristic is active
         *
         */
        bool is_lag_behind_active() {
            return _lag_behind;
        }

        /**
         * @brief Return whether the density lag heuristic is active
         *
         */
        bool is_lag_density_active() {
            return _lag_density;
        }

        /**
         * @brief Return whether the advancement density is active
         *
         */
        bool is_advancement_density_active() {
            return _advancement_density;
        }

        /**
         * @brief Return maximum number of steps
         *
         */
        int max_steps() {
            return _max_steps;
        }

    private:
        // Used heuristics
        bool _steps_limit;
        bool _lag_behind;
        bool _lag_density;
        bool _advancement_density;

        // General data
        int _max_offset;

        // Lag-behind + density
        int _min_offset_increase;
        int _scores_to_affect; // Once separation is found, the heuristic affects
                               // scores in the range [score, score + score_scope]
        bool _advanced_enough;
        int _K;
        std::vector<int> _last_max_offsets; // Vector of the maximum offsets in
                                            // the last K = (lag_threshold*0.3)*(penalties.mism()) + 1

        // Lag-behind heuristic
        int _lag_threshold;        // Maximum offset difference between max and curr offset (>= 10)

        // Advancement density
        int    _s_min;             // Minimum score to test the advancement density heuristic (>= 5/100*seq_len)
        double _density_threshold; // Threshold to drop the alignment

        // Maximum number of Theseus steps (1 step = 1 score)
        int  _max_steps;
    };

} // namespace theseus