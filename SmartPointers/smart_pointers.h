#include <iostream>

template<typename T>
class WeakPtr;

struct BaseControlBlock{
    size_t sharedCount = 0;
    size_t weakCount = 0;
    virtual void destroyObject() = 0;
    virtual void deallocate() = 0;
    virtual void* get() = 0;
    virtual ~BaseControlBlock() = default;
};

template<typename U, typename Deleter, typename Allocator>
struct RegularControlBlock: BaseControlBlock{
    U* ptr;
    Deleter deleter;
    Allocator alloc;

    RegularControlBlock(U* p, Deleter d, Allocator a): ptr(p), deleter(d), alloc(a){};

    void destroyObject() override{
        deleter(ptr);
    }

    void deallocate() override{
        using BlockAlloc = typename std::allocator_traits<Allocator>::template rebind_alloc<RegularControlBlock>;
        BlockAlloc blockAlloc = alloc;
        std::allocator_traits<BlockAlloc>::deallocate(blockAlloc, this, 1);
    }

    void* get(){
        return static_cast<void*>(ptr);
    }
};

template<typename U,typename Allocator>
struct MakerControlBlock: BaseControlBlock{
    alignas(U) char object[sizeof(U)];
    Allocator alloc;

    MakerControlBlock(Allocator a): alloc(a){};

    void destroyObject() override{
        using ObjAlloc = typename std::allocator_traits<Allocator>::template rebind_alloc<U>;
        ObjAlloc objAlloc = alloc;
        std::allocator_traits<ObjAlloc>::destroy(objAlloc, reinterpret_cast<U*>(object));
    }

    void deallocate() override{
        using BlockAlloc = typename std::allocator_traits<Allocator>::template rebind_alloc<MakerControlBlock>;
        BlockAlloc blockAlloc = alloc;
        std::allocator_traits<BlockAlloc>::deallocate(blockAlloc, this, 1);
    }

    void* get(){
        return reinterpret_cast<void*>(object);
    }
};

template <typename T>
class SharedPtr{
private:

    template<typename Y, typename Allocator, typename... Args>
    friend SharedPtr<Y> allocateShared(const Allocator& allocator, Args&&... args);

    template<typename Y>
    friend class WeakPtr;

    template<typename Y>
    friend class SharedPtr;

    T* ptr = nullptr;
    BaseControlBlock* controlBlock = nullptr;

    template<typename Allocator, typename... Args >
    SharedPtr(const Allocator& alloc, Args&&... args ){
        using BlockAlloc = typename std::allocator_traits<Allocator>::template rebind_alloc<MakerControlBlock<T, Allocator>>;
        BlockAlloc blockAlloc = alloc;
        MakerControlBlock<T, Allocator>* allocated_ptr = std::allocator_traits<BlockAlloc>::allocate(blockAlloc, 1);
        try {
            std::allocator_traits<BlockAlloc>::construct(blockAlloc, allocated_ptr, alloc);
        }
        catch (...) {
            std::allocator_traits<BlockAlloc>::deallocate(blockAlloc, allocated_ptr, 1);
            throw;
        }

        new(reinterpret_cast<T*>(allocated_ptr->object)) T(std::forward<Args>(args)...);

        allocated_ptr->sharedCount = 1;

        controlBlock = allocated_ptr;
    }

    SharedPtr(const WeakPtr<T>& other): ptr(other.ptr), controlBlock(static_cast<BaseControlBlock*>(other.controlBlock)) {
        ++controlBlock->sharedCount;
    }
public:
    SharedPtr() = default;

    template <typename Y, typename Deleter, typename Allocator>
    SharedPtr(Y* objPtr, Deleter deleter, Allocator alloc){
        ptr = objPtr;

        using BlockAlloc = typename std::allocator_traits<Allocator>::template rebind_alloc<RegularControlBlock<T,Deleter,Allocator>>;
        BlockAlloc blockAlloc = alloc;
        RegularControlBlock<T, Deleter, Allocator>* allocated_ptr = std::allocator_traits<BlockAlloc>::allocate(blockAlloc, 1);
        new(allocated_ptr) RegularControlBlock<T, Deleter, Allocator>(objPtr, std::move(deleter), std::move(alloc));

        allocated_ptr->sharedCount = 1;

        controlBlock = allocated_ptr;
    }

    template <typename Y, typename Deleter>
    SharedPtr(Y* objPtr, Deleter deleter): SharedPtr(objPtr, deleter, std::allocator<Y>()) {};

    template <typename Y>
    explicit SharedPtr(Y* objPtr): SharedPtr(objPtr, std::default_delete<Y>(), std::allocator<Y>()) {};


    template<typename Y>
    SharedPtr(const SharedPtr<Y>& other): ptr(other.ptr), controlBlock(other.controlBlock) {
        if(controlBlock) ++controlBlock->sharedCount;
    }

    SharedPtr(const SharedPtr& other): ptr(other.ptr), controlBlock(other.controlBlock) {
        if(controlBlock) ++controlBlock->sharedCount;
    }

    template<typename Y>
    SharedPtr(SharedPtr<Y>&& other): ptr(other.ptr), controlBlock(other.controlBlock) {
        other.ptr = nullptr;
        other.controlBlock = nullptr;
    }

    SharedPtr(SharedPtr&& other): ptr(other.ptr), controlBlock(other.controlBlock) {
        other.ptr = nullptr;
        other.controlBlock = nullptr;
    }

    template<typename Y>
    SharedPtr& operator=(const SharedPtr<Y>& other) {
        SharedPtr copy = other;
        swap(copy);
        return *this;
    }

    SharedPtr& operator=(const SharedPtr& other) {
        SharedPtr copy = other;
        swap(copy);
        return *this;
    }

    template<typename Y>
    SharedPtr& operator=(SharedPtr<Y>&& other) {
        SharedPtr copy = std::move(other);
        swap(copy);
        return *this;
    }

    SharedPtr& operator=(SharedPtr&& other) {
        SharedPtr copy = std::move(other);
        swap(copy);
        return *this;
    }

    size_t use_count() const{
        return controlBlock->sharedCount;
    }

    void swap(SharedPtr& other) {
        std::swap(ptr, other.ptr);
        std::swap(controlBlock, other.controlBlock);
    }

    void reset() {
        SharedPtr<T> other;
        swap(other);
    }

    template<typename Y>
    void reset(Y* object) {
        SharedPtr<T> other(object);
        swap(other);
    }

    T& operator*() const {
        if(ptr != nullptr) return *ptr;
        return *static_cast<T*>(controlBlock->get());
    }

    T* operator->() const {
        if(controlBlock == nullptr) return nullptr;
        if(ptr != nullptr) return ptr;
        return static_cast<T*>(controlBlock->get());
    }

    T* get() const {
        if(controlBlock == nullptr) return nullptr;
        if(ptr != nullptr) return ptr;
        return static_cast<T*>(controlBlock->get());
    }

    ~SharedPtr() {
        if(!controlBlock) return;
        --controlBlock->sharedCount;
        if(controlBlock->sharedCount > 0) return;
        controlBlock->destroyObject();
        if(controlBlock->weakCount == 0) controlBlock->deallocate();
    }
};


template<typename T, typename Allocator, typename... Args >
SharedPtr<T> allocateShared(const Allocator& alloc, Args&&... args ){
    return SharedPtr<T>(alloc, std::forward<Args>(args)...);
}

template<typename T, typename... Args >
SharedPtr<T> makeShared(Args&&... args ){
    return allocateShared<T, std::allocator<T>, Args...>(std::allocator<T>(), std::forward<Args>(args)...);
}


template <typename T>
class WeakPtr{
private:
    template<typename Y>
    friend class SharedPtr;

    template<typename Y>
    friend class WeakPtr;

    T* ptr = nullptr;
    BaseControlBlock* controlBlock = nullptr;
public:
    WeakPtr() = default;

    template<typename Y>
    WeakPtr(const SharedPtr<Y>& other): ptr(other.ptr), controlBlock(other.controlBlock) {
        ++controlBlock->weakCount;
    }


    template<typename Y>
    WeakPtr(const WeakPtr<Y>& other): ptr(other.ptr), controlBlock(other.controlBlock) {
        ++controlBlock->weakCount;
    }

    WeakPtr(const WeakPtr& other): ptr(other.ptr), controlBlock(other.controlBlock) {
        ++controlBlock->weakCount;
    }

    template<typename Y>
    WeakPtr(WeakPtr<Y>&& other): ptr(other.ptr), controlBlock(other.controlBlock) {
        other.ptr = nullptr;
        other.controlBlock = nullptr;
    }

    WeakPtr(WeakPtr&& other): ptr(other.ptr), controlBlock(other.controlBlock) {
        other.ptr = nullptr;
        other.controlBlock = nullptr;
    }

    template<typename Y>
    WeakPtr& operator=(const WeakPtr<Y>& other) {
        WeakPtr copy = other;
        swap(copy);
        return *this;
    }


    WeakPtr& operator=(const WeakPtr& other) {
        WeakPtr copy = other;
        swap(copy);
        return *this;
    }

    template<typename Y>
    WeakPtr& operator=(WeakPtr<Y>&& other) {
        WeakPtr copy = std::move(other);
        swap(copy);
        return *this;
    }

    WeakPtr& operator=(WeakPtr&& other) {
        WeakPtr copy = std::move(other);
        swap(copy);
        return *this;
    }

    size_t use_count() const{
        return controlBlock->sharedCount;
    }

    void swap(WeakPtr& other) {
        std::swap(ptr, other.ptr);
        std::swap(controlBlock, other.controlBlock);
    }

    bool expired() const {
        return controlBlock->sharedCount == 0;
    }

    SharedPtr<T> lock() const {
        if (!expired()) {
            return SharedPtr<T>(*this);
        }
        return SharedPtr<T>();
    }


    ~WeakPtr() {
        if(!controlBlock) return;
        --controlBlock->weakCount;
        if(controlBlock->weakCount > 0) return;
        if(controlBlock->sharedCount == 0) controlBlock->deallocate();
    }
};
