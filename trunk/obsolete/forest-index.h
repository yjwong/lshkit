/* 
    Copyright (C) 2008 Wei Dong <wdong@princeton.edu>. All Rights Reserved.
  
    This file is part of LSHKIT.
  
    LSHKIT is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LSHKIT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with LSHKIT.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef __LSHKIT_FLAT__
#define __LSHKIT_FLAT__

#include <algorithm>
#include <lshkit/topk.h>

namespace lshkit {

/// LSH Forest under development.
template <typename LSH, typename ACCESSOR, typename METRIC>
class ForestIndex
{
    BOOST_CONCEPT_ASSERT((LshConcept<LSH>));
public:
    typedef typename LSH::Parameter Parameter;
    typedef typename LSH::Domain Domain;
    typedef typename ACCESSOR::Key Key;

private:

    class Tree
    {
        std::vector<LSH> lsh_;
        unsigned min_, max_;

        class Node
        {
            
            ForestIndex *forest_;
            Tree *tree_;
            unsigned depth_;
            Node *parent_;
            std::vector<Node *> children_;
            std::vector<Key> data_;

        public:
            Node (ForestIndex *forest, Tree *tree, unsigned depth, Node *parent = 0)
                : forest_(forest), tree_(tree), depth_(depth), parent_(parent)
            {
            };

            ~Node ()
            {
                for (typename std::vector<Node *>::iterator it = children_.begin(); it != children_.end(); ++it)
                    if (*it != 0) delete *it;
            }

            bool empty ()
            {
                return data_.size() == 0 && children_.size() == 0;
            }

            void insert (Key key)
            {
                if (children_.size > 0)
                {
                    unsigned index = tree_->lsh_[depth_](forest_->accessor_(key));
                    children_[index]->insert(key);
                }
                else 
                {
                    data_.push_back(key);
                    if (depth_ < tree_->lsh_.size() && data_.size >= tree_->max_)
                    {
                        assert(children_.size() == 0);
                        LSH &lsh = tree_->lsh[depth_];
                        children_.resize(lsh.getRange());
                        for (typename std::vector<Node *>::iterator it = children_.begin(); it != children_.end(); ++it)
                        {
                            *it = new Node(forest_, tree_, depth_ + 1, this);
                        }
                        for (typename std::vector<Key>::iterator it = data_.begin(); it != data_.end(); ++it)
                        {
                            unsigned index = lsh(forest_->accessor_(*it));
                            children_[index]->insert(*it);
                        }
                        data_.clear();
                    }
                }
            }

            /*
            void prune ()
            {
                for (std::vector<Node *>::iterator it = children.begin(); it != children.end(); ++it)
                {
                    if (*it != 0) 
                    {
                        (*it)->prune();
                        if ((*it)->empty)
                        {
                            delete *it;
                            *it = 0;
                        }
                    }
                }

                if (children_->size() == 0)
                {
                    if (data_
                }
            }
            */

        } *root_;

        public:

        Tree ()
        {
        }

        template <typename ENGINE>
        void reset (const Parameter &param, ENGINE &engine, ForestIndex *forest, unsigned D, unsigned min, unsigned max)
        {
            lsh_.resize(D);
            for (unsigned i = 0; i < D; ++i)
            {
                lsh_[i].reset(param, engine);
            }
            root_ = new Node(forest, this, 0);
        }

        ~Tree ()
        {
            if (root_ != 0) delete root_;
        }

        void insert (Key key)
        {
            root_->insert(key);
        }
    };

    std::vector<Tree> trees_;

    ACCESSOR accessor_;
    METRIC metric_;
    

public:
    template <typename Engine>
    ForestIndex(const Parameter &param, Engine &engine, ACCESSOR accessor, METRIC metric, unsigned L, unsigned D, unsigned min, unsigned max) : 
        trees_(L), accessor_(accessor), metric_(metric)
    {
        for (unsigned i = 0; i < L; ++i)
        {
            trees_[i].reset(param, engine, this, D, min, max);
        }
    }

    void insert (Key key)
    {
        for (unsigned i = 0; i < trees_.size(); ++i)
        {
            trees_[i].insert(key);
        }
    }

    void query (const Domain &obj, Topk<Key> &topk, unsigned L = 0)
    {
        assert(0);
    }
};


}

#endif

