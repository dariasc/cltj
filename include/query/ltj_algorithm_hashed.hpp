/*
 * ltj_algorithm_hashed.hpp
 * Copyright (C) 2025 Author Diego Arias
 *
 *
 * This is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */



#ifndef RING_LTJ_ALGORITHM_HASHED_HPP
#define RING_LTJ_ALGORITHM_HASHED_HPP


#include <triple_pattern.hpp>
#include <dict/dict_map.hpp>
#include <query/ltj_iterator_basic.hpp>
#include <query/ltj_iterator_lite.hpp>
#include <veo/veo_simple.hpp>
#include <veo/veo_adaptive.hpp>
#include <results/results.hpp>
#include <util/rdf_util.hpp>

#include "ltj_iterator_hashed.hpp"
#include "hashed/cltj_hashed_index.h"

#define EXPT_TIME_SOL 0

namespace ltj {

    template<class iterator_t = ltj_iterator_hashed<cltj::hashed_ltj, uint8_t, uint64_t>,
             class veo_t = veo::veo_adaptive<iterator_t, util::trait_size> >
    class ltj_algorithm_hashed {

    public:
        typedef uint64_t value_type;
        typedef uint64_t size_type;
        typedef iterator_t ltj_iter_type;
        typedef typename ltj_iter_type::var_type var_type;
        typedef typename ltj_iter_type::index_scheme_type index_scheme_type;
        typedef typename ltj_iter_type::value_type const_type;
        typedef veo_t veo_type;
        typedef unordered_map<var_type, vector<ltj_iter_type*>> var_to_iterators_type;
        typedef vector<pair<var_type, const_type>> tuple_type;
        typedef vector<std::string> tuple_str_type;
        typedef chrono::high_resolution_clock::time_point time_point_type;
        //typedef ::util::results_collector<tuple_type> results_type;

    private:
        const vector<triple_pattern>* m_ptr_triple_patterns;
        veo_type m_veo;
        index_scheme_type* m_ptr_ring;
        vector<ltj_iter_type> m_iterators;
        var_to_iterators_type m_var_to_iterators;
        bool m_is_empty = false;


        void copy(const ltj_algorithm_hashed &o) {
            m_ptr_triple_patterns = o.m_ptr_triple_patterns;
            m_veo = o.m_veo;
            m_ptr_ring = o.m_ptr_ring;
            m_iterators = o.m_iterators;
            m_var_to_iterators = o.m_var_to_iterators;
            m_is_empty = o.m_is_empty;
        }


        inline void add_var_to_iterator(const var_type var, ltj_iter_type* ptr_iterator){
            auto it =  m_var_to_iterators.find(var);
            if(it != m_var_to_iterators.end()){
                it->second.push_back(ptr_iterator);
            }else{
                vector<ltj_iter_type*> vec = {ptr_iterator};
                m_var_to_iterators.insert({var, vec});
            }
        }

        void from_id_to_str(tuple_type &t, tuple_str_type &t_str, const std::vector<bool> &in_p,
                         dict::basic_map &dict_so, dict::basic_map &dict_p) {
              for(auto &p : t) {
                  if(in_p[p.first]) {
                      t_str[p.first] = dict_p.extract(p.second);
                  }else {
                      t_str[p.first] = dict_so.extract(p.second);
                  }
              }
        }

    public:


        ltj_algorithm_hashed() = default;

        ltj_algorithm_hashed(const vector<triple_pattern>* triple_patterns, index_scheme_type* ring){

            m_ptr_triple_patterns = triple_patterns;
            m_ptr_ring = ring;

            size_type i = 0;
            m_iterators.reserve(m_ptr_triple_patterns->size());
            for(const auto& triple : *m_ptr_triple_patterns){
                //Bulding iterators
                m_iterators.emplace_back(ltj_iter_type(&triple, m_ptr_ring));
                if(m_iterators[i].is_empty()){
                    m_is_empty = true;
                    return;
                }

                //For each variable we add the pointers to its iterators
                if(triple.o_is_variable()){
                    add_var_to_iterator(triple.term_o.value, &(m_iterators[i]));
                }
                if(triple.p_is_variable()){
                    add_var_to_iterator(triple.term_p.value, &(m_iterators[i]));
                }
                if(triple.s_is_variable()){
                    add_var_to_iterator(triple.term_s.value, &(m_iterators[i]));
                }
                ++i;
            }

            m_veo = veo_type(m_ptr_triple_patterns, &m_iterators, &m_var_to_iterators, m_ptr_ring);

        }

        //! Copy constructor
        ltj_algorithm_hashed(const ltj_algorithm_hashed &o) {
            copy(o);
        }

        //! Move constructor
        ltj_algorithm_hashed(ltj_algorithm_hashed &&o) {
            *this = move(o);
        }

        //! Copy Operator=
        ltj_algorithm_hashed &operator=(const ltj_algorithm_hashed &o) {
            if (this != &o) {
                copy(o);
            }
            return *this;
        }

        //! Move Operator=
        ltj_algorithm_hashed &operator=(ltj_algorithm_hashed &&o) {
            if (this != &o) {
                m_ptr_triple_patterns = move(o.m_ptr_triple_patterns);
                m_veo = move(o.m_veo);
                m_ptr_ring = move(o.m_ptr_ring);
                m_iterators = move(o.m_iterators);
                m_var_to_iterators = move(o.m_var_to_iterators);
                m_is_empty = o.m_is_empty;
            }
            return *this;
        }

        void swap(ltj_algorithm_hashed &o) {
            std::swap(m_ptr_triple_patterns, o.m_ptr_triple_patterns);
            std::swap(m_veo, o.m_veo);
            std::swap(m_ptr_ring, o.m_ptr_ring);
            std::swap(m_iterators, o.m_iterators);
            std::swap(m_var_to_iterators, o.m_var_to_iterators);
            std::swap(m_is_empty, o.m_is_empty);
        }


        /**
        *
        * @param res               Results
        * @param limit_results     Limit of results
        * @param timeout_seconds   Timeout in seconds
        */
        template<class results_type>
        void join(/*vector<tuple_type> &res,*/
                  results_type &res,
                  const size_type limit_results = 0, const size_type timeout_seconds = 0){
            if(m_is_empty) return;
            time_point_type start = chrono::high_resolution_clock::now();
            tuple_type t(m_veo.size());
            search(0, t, res, start, limit_results, timeout_seconds);
        };

        /**
        *
        * @param res               Results
        * @param limit_results     Limit of results
        * @param timeout_seconds   Timeout in seconds
        */
        template<class results_type>
        void join_str(/*vector<tuple_type> &res,*/
                  results_type &res,
                  const std::vector<bool> &in_p,
                  dict::basic_map &dict_so, dict::basic_map &dict_p,
                  const size_type limit_results = 0, const size_type timeout_seconds = 0){
        };

        /**
         *
         * @param j                 Index of the variable
         * @param tuple             Tuple of the current search
         * @param res               Results
         * @param start             Initial time to check timeout
         * @param limit_results     Limit of results
         * @param timeout_seconds   Timeout in seconds
         */
        template<class results_type>
        bool search_str(const size_type j, tuple_type &tuple, tuple_str_type &t_str,
                        results_type &res,
                        const std::vector<bool> &in_p,
                        dict::basic_map &dict_so, dict::basic_map &dict_p,
                        const time_point_type start,
                        const size_type limit_results = 0, const size_type timeout_seconds = 0){
        };

        /**
         *
         * @param j                 Index of the variable
         * @param tuple             Tuple of the current search
         * @param res               Results
         * @param start             Initial time to check timeout
         * @param limit_results     Limit of results
         * @param timeout_seconds   Timeout in seconds
         */
        template<class results_type>
        bool search(const size_type j, tuple_type &tuple, results_type &res,
                    const time_point_type start,
                    const size_type limit_results = 0, const size_type timeout_seconds = 0){
            //(Optional) Check timeout
            if(timeout_seconds > 0){
                time_point_type stop = chrono::high_resolution_clock::now();
                auto sec = chrono::duration_cast<chrono::seconds>(stop-start).count();
                if(sec > timeout_seconds) {
                    return false;
                }
            }

            //(Optional) Check limit
            if(limit_results > 0 && res.size() == limit_results) {
                return false;
            }

            if(j == m_veo.size()){
                //Report results
                //res.emplace_back(tuple);
                res.add(tuple);
#if EXPT_TIME_SOL
                if(res.size() % 1000 == 0){
                    time_point_type stop = chrono::high_resolution_clock::now();
                    auto sec = chrono::duration_cast<chrono::seconds>(stop-start).count();
                    std::cerr << res.size() << ";" << sec << std::endl;
                }
#endif
            }else{
                var_type x_j = m_veo.next();
                //cout << "Variable: " << (uint64_t) x_j << endl;
                //cout << (uint64_t) x_j << endl;
                vector<ltj_iter_type*>& itrs = m_var_to_iterators[x_j];
                bool ok;
                if(itrs.size() == 1 && itrs[0]->in_last_level()) {//Lonely variables
                    //cout << "Seeking (last level)" << endl;
                    auto cnt = itrs[0]->children_var(x_j);

                    //cout << "Results: " << results.size() << endl;
                    //cout << "Seek (last level): (" << (uint64_t) x_j << ": size=" << results.size() << ")" <<endl;
                    for (size_type k = 1; k <= cnt; ++k) {
                        auto c = itrs[0]->seek(x_j, k);
                        //1. Adding result to tuple
                        tuple[j] = {x_j, c};
                        //2. Going down in the trie by setting x_j = c (\mu(t_i) in paper)
                        itrs[0]->down_var(x_j, c);
                        m_veo.down();
                        //2. Search with the next variable x_{j+1}
                        ok = search(j + 1, tuple, res, start, limit_results, timeout_seconds);
                        if(!ok) return false;
                        //4. Going up in the trie by removing x_j = c
                        itrs[0]->up();
                        m_veo.up();
                    }
                }else {
                    sort(itrs.begin(), itrs.end(), [&](auto &a, auto &b) {
                        return a->children_var(x_j) < b->children_var(x_j);
                    });
                    auto cnt = itrs[0]->children_var(x_j);
                    for (size_type k = 1; k <= cnt; ++k) {
                        auto c = itrs[0]->seek(x_j, k);
                        bool matched = true;
                        for (size_type i = 1; i < itrs.size(); ++i) {
                            if (!itrs[i]->exists_var(x_j, c)) {
                                matched = false;
                                break;
                            }
                        }
                        if (matched) {
                            //1. Adding result to tuple
                            tuple[j] = {x_j, c};
                            //2. Going down in the trie by setting x_j = c (\mu(t_i) in paper)
                            for (ltj_iter_type* iter : itrs) {
                                iter->down_var(x_j, c);
                            }
                            m_veo.down();
                            //2. Search with the next variable x_{j+1}
                            ok = search(j + 1, tuple, res, start, limit_results, timeout_seconds);
                            if(!ok) return false;
                            //4. Going up in the trie by removing x_j = c
                            for (ltj_iter_type* iter : itrs) {
                                iter->up();
                            }
                            m_veo.up();
                        }
                    }
                }
                m_veo.done();
            }
            return true;
        };

        void print_veo(unordered_map<uint8_t, string> &ht){
            cout << "veo: ";
            for(uint64_t j = 0; j < m_veo.size(); ++j){
                cout << "?" << ht[m_veo.next()] << " ";
            }
            cout << endl;
        }

        void print_query(unordered_map<uint8_t, string> &ht){
            cout << "Query: " << endl;
            for(size_type i = 0; i <  m_ptr_triple_patterns->size(); ++i){
                m_ptr_triple_patterns->at(i).print(ht);
                if(i < m_ptr_triple_patterns->size()-1){
                    cout << " . ";
                }
            }
            cout << endl;
        }

        void print_results(vector<tuple_type> &res, unordered_map<uint8_t, string> &ht){
            cout << "Results: " << endl;
            uint64_t i = 1;
            for(tuple_type &tuple :  res){
                cout << "[" << i << "]: ";
                for(pair<var_type, value_type> &pair : tuple){
                    cout << "?" << ht[pair.first] << "=" << pair.second << " ";
                }
                cout << endl;
                ++i;
            }
        }


    };
}

#endif // RING_LTJ_ALGORITHM_HASHED_HPP