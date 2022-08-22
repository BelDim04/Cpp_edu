#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>


template<typename T, typename Allocator = std::allocator<T>>
class List {
public:
    template<typename U, typename Alloc>
    friend bool operator==(const List<U, Alloc> &a, const List<U, Alloc> &b);

    struct BaseNode {
        BaseNode *next = this;
        BaseNode *prev = this;
    };

    struct Node : public BaseNode {
        T value;

        Node() = default;

        Node(const T &val) : value(val) {}
    };

    template<bool IsConst>
    struct common_iterator {
        BaseNode* ptr;
        using pointer = std::conditional_t<IsConst, const T *, T *>;
        using reference = std::conditional_t<IsConst, const T &, T &>;
        using value_type = std::conditional_t<IsConst, const T, T>;
        using difference_type = std::ptrdiff_t;
        using iterator_category = std::bidirectional_iterator_tag;

        common_iterator<IsConst>(BaseNode *v) : ptr(v) {};

        common_iterator<IsConst>(const common_iterator<false> &it) : ptr(it.ptr) {};

        common_iterator<IsConst> &operator=(const common_iterator<false> &it) {
            ptr = it.ptr;
            return *this;
        };

        common_iterator<IsConst> &operator--() {
            ptr = ptr->prev;
            return *this;
        };

        common_iterator<IsConst> &operator++() {
            ptr = ptr->next;
            return *this;
        };

        const common_iterator<IsConst> operator++(int) {
            common_iterator <IsConst> temp = *this;
            ptr = ptr->next;
            return temp;
        };

        const common_iterator<IsConst> operator--(int) {
            common_iterator <IsConst> temp = *this;
            ptr = ptr->prev;
            return temp;
        };

        std::conditional_t<IsConst, const T &, T &> operator*() const{
            return static_cast<Node *>(ptr)->value;
        };

        template<bool IsOtherConst>
        bool operator==(const common_iterator<IsOtherConst> other) {
            return ptr == other.ptr;
        }

        template<bool IsOtherConst>
        bool operator!=(const common_iterator<IsOtherConst> other) {
            return !(*this == other);
        }
    };

    Allocator alloc;
    using NodeAlloc = typename std::allocator_traits<Allocator>::template rebind_alloc<Node>;
    NodeAlloc nodeAlloc;

    BaseNode fakeNodeObject = BaseNode();
    BaseNode *fakeNode = &fakeNodeObject;
    size_t sz = 0;

    using value_type = T;
    using reference = T &;
    using const_reference = const T &;
    using size_type = size_t;
    using const_iterator = common_iterator<true>;
    using iterator = common_iterator<false>;
    using const_reverse_iterator = std::reverse_iterator<common_iterator<true>>;
    using reverse_iterator = std::reverse_iterator<common_iterator<false>>;
    using allocator_type = Allocator;


    List<T, Allocator>() = default;

    template <typename ...Args>
    List<T, Allocator>(size_t n, Args ...args, Allocator a) : alloc(a), nodeAlloc(alloc) {
        for (size_t i = 0; i < n; ++i) {
            try {
                this->emplace(end(), args...);
            }
            catch (...) {
                for (size_t j = 0; j < i; ++j) {
                    this->erase(--end());
                }
                throw;
            }
        }
    };

    List<T, Allocator>(size_t n): List<T, Allocator>(n, alloc) {}

    List<T, Allocator>(size_t n, const T &val) : List<T, Allocator>(n, val, alloc){}

    List<T, Allocator>(Allocator a) : alloc(a), nodeAlloc(alloc) {};

    List<T, Allocator>(const List<T, Allocator> &other) : alloc(
            std::allocator_traits<Allocator>::select_on_container_copy_construction(other.alloc)), nodeAlloc(
            std::allocator_traits<Allocator>::select_on_container_copy_construction(other.nodeAlloc)) {
        iterator it = other.begin();
        while (sz != other.sz) {
            try {
                this->emplace(end(), *it);
                ++it;
            }
            catch (...) {
                while (sz != 0) {
                    this->erase(--end());
                }
                throw;
            }
        }
    };

    List<T, Allocator> &operator=(const List<T, Allocator> &other) {
        Allocator alloc_copy = alloc;
        NodeAlloc nodeAlloc_copy = nodeAlloc;
        if (std::allocator_traits<Allocator>::propagate_on_container_copy_assignment::value) {
            alloc = other.alloc;
            nodeAlloc = other.nodeAlloc;
        }
        size_t beforeAssignmentSz = sz;
        iterator it = other.begin();
        while (sz - beforeAssignmentSz != other.sz) {
            try {
                this->emplace(end(), *it);
                ++it;
            }
            catch (...) {
                while (sz > beforeAssignmentSz) {
                    this->erase(--end());
                }
                throw;
            }
        }
        alloc = alloc_copy;
        nodeAlloc = nodeAlloc_copy;
        while (sz > other.sz) {
            this->erase(begin());
        }
        if (std::allocator_traits<Allocator>::propagate_on_container_copy_assignment::value) {
            alloc = other.alloc;
            nodeAlloc = other.nodeAlloc;
        }
        return *this;
    };


    List(List<T, Allocator>&& other) {
        alloc = std::move(other.alloc);
        nodeAlloc = std::move(other.nodeAlloc);
        fakeNodeObject = std::move(other.fakeNodeObject);
        fakeNode = &fakeNodeObject;
        fakeNode->next->prev = fakeNode;
        fakeNode->prev->next = fakeNode;
        other.fakeNode = nullptr;
        sz = other.sz;
        other.sz = 0;
    }

    List& operator=(List<T, Allocator>&& other) {
        alloc = std::move(other.alloc);
        nodeAlloc = std::move(other.nodeAlloc);
        fakeNodeObject = std::move(other.fakeNodeObject);
        fakeNode = &fakeNodeObject;
        fakeNode->next->prev = fakeNode;
        fakeNode->prev->next = fakeNode;
        other.fakeNode = nullptr;
        sz = other.sz;
        other.sz = 0;
        return *this;
    }

    ~List<T, Allocator>() {
        while (sz > 0) {
            this->pop_back();
        }
    };


    iterator begin() const {
        return iterator(fakeNode->next);
    }

    const_iterator cbegin() const {
        return const_iterator(fakeNode->next);
    }

    reverse_iterator rbegin() const {
        return reverse_iterator(fakeNode);
    }

    const_reverse_iterator crbegin() const {
        return const_reverse_iterator(fakeNode);
    }

    iterator end() const {
        return iterator(fakeNode);
    }

    const_iterator cend() const {
        return const_iterator(fakeNode);
    }

    reverse_iterator rend() const {
        return reverse_iterator(fakeNode->next);
    }

    const_reverse_iterator crend() const {
        return const_reverse_iterator(fakeNode->next);
    }

    Allocator get_allocator() const {
        return alloc;
    }

    size_t size() const {
        return sz;
    }

    void insert(List<T, Allocator>::common_iterator<true> pos, const T &value) {
        try{
            emplace(pos, value);
        }
        catch (...){
            throw;
        }
    };

    void insert(List<T, Allocator>::common_iterator<true> pos, T&& value) {
        try{
            emplace(pos, std::move(value));
        }
        catch (...){
            throw;
        }
    }

    template<typename...Args>
    iterator emplace(List<T, Allocator>::common_iterator<true> pos, Args&& ...args) {
        Node *allocated_ptr = std::allocator_traits<NodeAlloc>::allocate(nodeAlloc, 1);
        try {
            std::allocator_traits<NodeAlloc>::construct(nodeAlloc, allocated_ptr, std::forward<Args>(args)...);
        }
        catch (...) {
            std::allocator_traits<NodeAlloc>::deallocate(nodeAlloc, allocated_ptr, 1);
            throw;
        }
        allocated_ptr->prev = pos.ptr->prev;
        allocated_ptr->next = pos.ptr;
        pos.ptr->prev->next = allocated_ptr;
        pos.ptr->prev = allocated_ptr;
        ++sz;
        return iterator(allocated_ptr);
    }

    void erase(List<T, Allocator>::common_iterator<true> pos) {
        pos.ptr->prev->next = pos.ptr->next;
        pos.ptr->next->prev = pos.ptr->prev;
        std::allocator_traits<NodeAlloc>::destroy(nodeAlloc, static_cast<Node *>(pos.ptr));
        std::allocator_traits<NodeAlloc>::deallocate(nodeAlloc, static_cast<Node *>(pos.ptr), 1);
        --sz;
    };

    void push_back(const T &value) {
        try {
            emplace(end(), value);
        }
        catch (...) {
            throw;
        }
    }

    void push_front(const T &value) {
        try {
            emplace(begin(), value);
        }
        catch (...) {
            throw;
        }
    }

    void pop_back() {
        erase(--end());
    }

    void pop_front() {
        erase(begin());
    }
};

template<typename T, typename Allocator>
bool operator==(const List<T, Allocator> &a, const List<T, Allocator> &b) {
    if (a.sz != b.sz) return false;
    typename List<T, Allocator>::iterator it1 = a.begin();
    typename List<T, Allocator>::iterator it2 = b.begin();
    while (it1 != a.end()) {
        if (*it1 != *it2) return false;
        ++it1;
        ++it2;
    }
    return true;
}

template<typename T, typename Allocator>
bool operator!=(const List<T, Allocator> &a, const List<T, Allocator> &b) {
    return !(a == b);
}




template<typename Key,typename Value,typename Hash=std::hash<Key>,typename Equal=std::equal_to<Key>,typename Alloc=std::allocator<std::pair<const Key, Value>>>
class UnorderedMap{
public:
    using NodeType = std::pair<const Key, Value>;

private:

    Alloc alloc;

    Hash keyHash;

    Equal keyEqual;

    struct ListNode{
        NodeType* element;
        size_t cached;

        ListNode(NodeType* e, size_t c): element(e), cached(c){}
    };
    using ListNodeAlloc = typename std::allocator_traits<Alloc>::template rebind_alloc<ListNode>;

    List<ListNode, ListNodeAlloc> elementsList;
    std::vector<typename List<ListNode, ListNodeAlloc>::iterator> keyBuckets = std::vector<typename List<ListNode, ListNodeAlloc>::iterator>(10, elementsList.end());

    template<bool IsConst>
    struct common_iterator {
        typename List<ListNode, ListNodeAlloc>::iterator listItr;
        using pointer = std::conditional_t<IsConst, const NodeType *, NodeType *>;
        using reference = std::conditional_t<IsConst, const NodeType &, NodeType &>;
        using value_type = std::conditional_t<IsConst, const NodeType, NodeType>;
        using difference_type = std::ptrdiff_t;
        using iterator_category = std::forward_iterator_tag;

        common_iterator<IsConst>(typename List<ListNode, ListNodeAlloc>::iterator itr) : listItr(itr) {};

        common_iterator<IsConst>(const common_iterator<false> &it) : listItr(it.listItr) {};

        common_iterator<IsConst> &operator=(const common_iterator<false> &it) {
            listItr = it.listItr;
            return *this;
        };

        common_iterator<IsConst> &operator++() {
            ++listItr;
            return *this;
        };

        const common_iterator<IsConst> operator++(int) {
            common_iterator <IsConst> temp = *this;
            ++listItr;
            return temp;
        };

        std::conditional_t<IsConst, const NodeType &, NodeType &> operator*() const{
            return *(*listItr).element;
        };

        std::conditional_t<IsConst, const NodeType*, NodeType*> operator->() {
            return (*listItr).element;
        }

        template<bool IsOtherConst>
        bool operator==(const common_iterator<IsOtherConst> other) {
            return listItr == other.listItr;
        }

        template<bool IsOtherConst>
        bool operator!=(const common_iterator<IsOtherConst> other) {
            return !(*this == other);
        }
    };

    float m_load_factor = 1.0;
public:

   using iterator = common_iterator<false>;
   using const_iterator = common_iterator<true>;


    UnorderedMap() = default;

    UnorderedMap(const UnorderedMap& other) {
        for (auto& x : other.elementsList) {
            emplace(*(x.element));
        }
    }

    UnorderedMap& operator=(const UnorderedMap& other) {
        UnorderedMap copy = *other;
        std::swap(*this, copy);
        return *this;
    }

    UnorderedMap& operator=(UnorderedMap&& other) {
        keyBuckets = std::move(other.keyBuckets);
        elementsList = std::move(other.elementsList);
        return *this;
    }

    UnorderedMap(UnorderedMap&& other) {
        keyBuckets = std::move(other.keyBuckets);
        elementsList = std::move(other.elementsList);
    }

    ~UnorderedMap(){
        while(elementsList.sz > 0) erase(begin());
    }

   size_t size() const{
       return elementsList.size();
   }

    iterator begin() const {
        return iterator(elementsList.begin());
    }

    const_iterator cbegin() const {
        return const_iterator(elementsList.begin());
    }

    iterator end() const {
        return iterator(elementsList.end());
    }

    const_iterator cend() const {
        return const_iterator(elementsList.end());
    }

    float max_load_factor() const{
        return m_load_factor;
    }

    void max_load_factor( float z ){
       m_load_factor = z;
   }

   float load_factor() const{
       return static_cast<float>(size()) / keyBuckets.size();
   }

   size_t max_size() const{
       return keyBuckets.size() * m_load_factor;
   }

   template<typename... Args>
   std::pair<iterator, bool> emplace(Args&&... args){
       if(size() + 1 > m_load_factor * keyBuckets.size()) rehash(keyBuckets.size() * 2);
       NodeType *allocated_ptr = std::allocator_traits<Alloc>::allocate(alloc, 1);
       std::allocator_traits<Alloc>::construct(alloc, allocated_ptr, std::forward<Args>(args)...);
       size_t hash = keyHash((*allocated_ptr).first);
       size_t index = hash % keyBuckets.size();
       iterator fIt = find_by_hash((*allocated_ptr).first, hash);
       if(fIt != end()) {
           std::allocator_traits<Alloc>::destroy(alloc, allocated_ptr);
           std::allocator_traits<Alloc>::deallocate(alloc, allocated_ptr, 1);
           return {fIt, false};
       }
       keyBuckets[index] = elementsList.emplace(keyBuckets[index], ListNode(allocated_ptr, hash));
       return {keyBuckets[index], true};
   }

   iterator find(const Key& key){
       size_t hash = keyHash(key);
       return find_by_hash(key, hash);
   }

   iterator find_by_hash(const Key& key, size_t hash){
       typename List<ListNode, ListNodeAlloc>::iterator it = keyBuckets[hash % keyBuckets.size()];
       while(it != elementsList.end() && (*it).cached % keyBuckets.size() == hash % keyBuckets.size()){
           if(keyEqual((*(*it).element).first, key)){
               return it;
           }
           ++it;
       }
       return end();
   }

   void rehash(size_t n){
       std::vector<typename List<ListNode, ListNodeAlloc>::iterator> nBuckets(n, elementsList.end());
       typename List<ListNode, ListNodeAlloc>::iterator it = elementsList.begin();
       typename List<ListNode, ListNodeAlloc>::BaseNode* curEnd = elementsList.fakeNode;
       while(it != elementsList.end()){
           size_t ind = (*it).cached % n;
           typename List<ListNode, ListNodeAlloc>::BaseNode* ptr = it.ptr;
           ++it;
           if(nBuckets[ind] == elementsList.end()){
               ptr->prev = curEnd;
               curEnd->next = ptr;
               ptr->next = elementsList.fakeNode;
               elementsList.fakeNode->prev = ptr;
               curEnd = ptr;
           }
           else{
               ptr->prev = nBuckets[ind].ptr->prev;
               nBuckets[ind].ptr->prev->next = ptr;
               ptr->next = nBuckets[ind].ptr;
               nBuckets[ind].ptr->prev = ptr;
           }
           nBuckets[ind] = typename List<ListNode, ListNodeAlloc>::iterator(ptr);
       }
       keyBuckets = nBuckets;
   }

    void erase(iterator it){
        typename List<ListNode, ListNodeAlloc>::iterator lIt = it.listItr;
        size_t pInd = (*lIt).cached % keyBuckets.size();
        if(keyBuckets[pInd] == lIt){
            typename List<ListNode, ListNodeAlloc>::iterator nextIt(lIt.ptr->next);
            if(nextIt != elementsList.end() && (*nextIt).cached % keyBuckets.size() == pInd) keyBuckets[pInd] = nextIt;
            else keyBuckets[pInd] = elementsList.end();
        }
        std::allocator_traits<Alloc>::destroy(alloc, (*lIt).element);
        std::allocator_traits<Alloc>::deallocate(alloc, (*lIt).element, 1);
        (*lIt).element = nullptr;
        elementsList.erase(lIt);
    }

    template<typename InputIterator>
    void erase(InputIterator itFirst, InputIterator itLast) {
        while (itFirst != itLast) {
            InputIterator nxt = itFirst;
            ++nxt;
            erase(itFirst);
            itFirst = nxt;
        }
    }

   void reserve(size_t n){
       if(n > m_load_factor * keyBuckets.size()) rehash(2 * n / m_load_factor);
   }

    std::pair<iterator, bool> insert(const NodeType& node) {
        return emplace(node);
    }

    std::pair<iterator, bool> insert(NodeType&& node) {
        return emplace(std::move(const_cast<Key&>(node.first)), std::move(node.second));
    }

    template<typename InputIterator>
    void insert(InputIterator itFirst, InputIterator itLast) {
        while (itFirst != itLast) {
            emplace(*itFirst);
            ++itFirst;
        }
    }

    Value& at(const Key& key){
       iterator fIt = find(key);
       if(fIt == end()) throw 1;
       else return (*fIt).second;
   }

    Value& operator[](const Key& key) {
        iterator fIt = find(key);
        if(fIt == end()) return emplace(std::pair(key, Value())).first->second;
        else return fIt->second;
    }

};
