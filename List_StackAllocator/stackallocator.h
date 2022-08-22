#include <iostream>
#include <list>
#include <memory>
#include <cstddef>
#include <type_traits>

template<size_t N>
class alignas(std::max_align_t) StackStorage {
private:
    char storage[N];
    size_t cur_position = 0;

public:
    StackStorage<N>() = default;

    StackStorage<N>(const StackStorage<N> &) = delete;

    StackStorage<N> operator=(const StackStorage<N> &) = delete;

    ~StackStorage<N>() = default;

    char *place_in_storage(size_t align, size_t num_of_bytes) {
        cur_position = (cur_position + align - 1) / align * align + num_of_bytes;
        return storage + cur_position - num_of_bytes;
    }
};

template<typename T, size_t N>
class StackAllocator {
    template<typename U, size_t M>
    friend bool operator==(const StackAllocator<U, M> &a, const StackAllocator<U, M> &b);

public:
    StackStorage<N> *stackStorage;

    using value_type = T;

    using propagate_on_container_copy_assignment = std::true_type;

    template<typename U>
    struct rebind {
        using other = StackAllocator<U, N>;
    };


    StackAllocator<T, N>() = delete;

    template<typename U>
    StackAllocator<T, N>(const StackAllocator<U, N> &other): stackStorage(other.stackStorage) {};

    StackAllocator<T, N>(StackStorage<N> &stackStorage) : stackStorage(&stackStorage) {};

    ~StackAllocator() = default;

    StackAllocator<T, N> &operator=(const StackAllocator<T, N> &other) {
        stackStorage = other.stackStorage;
        return *this;
    }

    StackAllocator<T, N> select_on_container_copy_construction() {
        return *this;
    }

    T *allocate(size_t n) {
        return reinterpret_cast<T *>(stackStorage->place_in_storage(alignof(T), sizeof(T) * n));
    }

    template<typename...Args>
    void construct(T *ptr, const Args &...args) {
        new(ptr) T(args...);
    }

    void destroy(T *ptr) {
        ptr->~T();
    }

    void deallocate(T *, size_t) {}
};

template<class T, size_t N>
bool operator==(const StackAllocator<T, N> &a, const StackAllocator<T, N> &b) {
    return a.stackStorage == b.stackStorage;
}

template<class T, size_t N>
bool operator!=(const StackAllocator<T, N> &a, const StackAllocator<T, N> &b) {
    return !(a == b);
}


template<typename T, typename Allocator = std::allocator<T>>
class List {
    template<typename U, typename Alloc>
    friend bool operator==(const List<U, Alloc> &a, const List<U, Alloc> &b);

private:
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

        common_iterator<IsConst> operator++(int) {
            common_iterator <IsConst> temp = *this;
            ptr = ptr->next;
            return temp;
        };

        common_iterator<IsConst> operator--(int) {
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

public:
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

    template<typename...Args>
    void emplace(List<T, Allocator>::common_iterator<true> pos, const Args &...args) {
        Node *allocated_ptr = std::allocator_traits<NodeAlloc>::allocate(nodeAlloc, 1);
        try {
            std::allocator_traits<NodeAlloc>::construct(nodeAlloc, allocated_ptr, args...);
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