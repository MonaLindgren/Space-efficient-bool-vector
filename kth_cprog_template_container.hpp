#include <cstdint>
#include <cstddef>
#include <initializer_list>
#include <stdexcept>
#include <type_traits>
#include <iostream>
#include <utility>
#include <limits>
#include <iterator>


//generell template:
namespace lab1{
template <typename T> 
class Vector
{
private:
    std::size_t _size;
    std::size_t _capacity;
    T * _front;
    void reallocate();
    std::size_t max(const std::size_t&, const std::size_t&) const;
public:
    template <typename U>
    friend void swap(Vector<U>&, Vector<U>&);
    Vector();
    explicit Vector(std::size_t);
    explicit Vector(std::size_t, T);
    Vector(std::initializer_list<T>);
    Vector(const Vector<T>&);
    Vector(Vector<T>&&);
    ~Vector();
    Vector& operator=(const Vector<T>&);
    Vector& operator=(Vector<T>&&);
    
    T& operator[](std::size_t); //indexing operator
    const T& operator[](std::size_t) const;
    void reset();
    std::size_t size() const;
    std::size_t capacity() const;
    void push_back(T);
    void insert(std::size_t, T);
    void erase(std::size_t);
    void clear();
    T* begin();
    const T* begin() const;
    T* end();
    const T* end()const;
    T* find(T const& elem);
    const T* find(T const& elem) const;
};
//specialisering vector<bool>
template <> 
class Vector<bool>
{
public:
    //måste forward-deklarera dessa för att kunna deklarera friend-operatorer
    class const_iterator;
    class iterator;
    friend class const_iterator;
    friend class iterator;
    friend bool operator==(const Vector<bool>::const_iterator&, const Vector<bool>::const_iterator&);
    friend bool operator!=(const Vector<bool>::const_iterator&, const Vector<bool>::const_iterator&);
    friend bool operator<(const Vector<bool>::const_iterator&, const Vector<bool>::const_iterator&);
    friend bool operator<=(const Vector<bool>::const_iterator&, const Vector<bool>::const_iterator&);
    friend bool operator>(const Vector<bool>::const_iterator&, const Vector<bool>::const_iterator&);
    friend bool operator>=(const Vector<bool>::const_iterator&, const Vector<bool>::const_iterator&);
    //friend long int operator -(const Vector<bool>::const_iterator&, const Vector<bool>::const_iterator&);
    
    class Reference {
        friend class Vector<bool>;
        Reference ();
        std::size_t _int_index; //vilken int det handlar om
        std::size_t _bit_index; //vilken bit i inten det handlar om...
        unsigned int * _front;
        public:
        ~Reference ();
        Reference& operator= (bool);
        Reference& operator=( const Reference& );
        operator bool () const;
    };
public:
    class const_iterator: public std::iterator<
                            std::random_access_iterator_tag,
                            Vector<bool>::Reference
                                        >{
        friend class Vector<bool>;//för att kunna anropa den privata konstruktorn.
        const Vector<bool>* _vec;
        
        protected:
            std::ptrdiff_t _index; //index in vector<bool> pointed to by iterator
        public:
            const_iterator(); 
            explicit const_iterator(const Vector<bool>&, const std::size_t&); 
            virtual ~const_iterator(){}
            const_iterator(const const_iterator&);
            const_iterator& operator=(const const_iterator&);
            virtual std::string type() const {return "const";}
            friend bool operator==(const const_iterator&, const const_iterator&);
            friend bool operator!=(const const_iterator&, const const_iterator&);
            friend bool operator<(const const_iterator&, const const_iterator&);
            friend bool operator<=(const const_iterator&, const const_iterator&);
            friend bool operator>(const const_iterator&, const const_iterator&);
            friend bool operator>=(const const_iterator&, const const_iterator&);
            friend std::ptrdiff_t operator-(const const_iterator&) const;
            const_iterator& operator++();
            const_iterator& operator+=(int);
            const_iterator& operator-=(int);
            const_iterator operator++(int);
            const_iterator& operator--();
            const_iterator operator--(int);
            const_iterator operator+(int);
            const_iterator operator-(int);
            const Vector<bool>::Reference operator*() const;
    };

    class iterator : public const_iterator {
        Vector<bool>* _vec;
        public:
            //using Vector<bool>::const_iterator::operator-;
        iterator();
        virtual ~iterator(){}
        explicit iterator(Vector<bool>&, const std::size_t&);
        iterator(const iterator&);
        iterator& operator=(const iterator&);
        virtual std::string type() const {return "nonconst";}
        iterator& operator++();
        iterator operator++(int);
        iterator& operator--();
        iterator operator--(int);
        iterator operator+(int);
        iterator operator-(int);
        Vector<bool>::Reference operator*();
    };

private:
    std::size_t _chunk = sizeof(unsigned int)*8;
    std::size_t _size;
    std::size_t _capacity;
    unsigned int * _front;
    void reallocate();
    std::size_t max(const std::size_t&, const std::size_t&) const;
    std::size_t size_convert(std::size_t) const;
    std::pair<std::size_t, std::size_t> intern_indexes(std::size_t) const;
public:
    friend void swap(Vector<bool>&, Vector<bool>&);
    unsigned int to_int() const;
	Vector();
    explicit Vector(std::size_t);
    explicit Vector(std::size_t, bool);
    Vector(std::initializer_list<bool>);
    Vector(const Vector<bool>&);
    Vector(Vector<bool>&&);
    ~Vector();
    Vector& operator=(const Vector<bool>&);
    Vector& operator=(Vector<bool>&&);
    
    Reference operator[](std::size_t);
    const Reference operator[](std::size_t) const;
    void reset();
    std::size_t size() const;
    std::size_t capacity() const;
    void push_back(bool);
    void clear();
    iterator begin();
    const_iterator begin() const;
    iterator end();
    const_iterator end() const;
    //T* find(T const& elem);
    //const T* find(T const& elem) const;
};

std::size_t Vector<bool>::size_convert(std::size_t x) const{
    std::size_t size = x/_chunk;
    size = (size*_chunk<x)?size+1:size;
    return size;
}

std::pair<std::size_t, std::size_t> Vector<bool>::intern_indexes(std::size_t index) const {
    return std::make_pair(index/_chunk, index%_chunk);
}

unsigned int Vector<bool>::to_int() const{
    if(_size>_chunk){
        throw std::out_of_range ("vector too big");
    } else {
        return _front[0];
    }
}
Vector<bool>::Reference::Reference(){
    _front = 0;
    _int_index = 0;
    _bit_index = 0;
}
Vector<bool>::Reference::~Reference(){
    _front = 0;
    _int_index = 0;
    _bit_index = 0;
}
Vector<bool>::Reference& Vector<bool>::Reference::operator= (bool x){
    if (_front!=0)
    {
        _front[_int_index] = (_front[_int_index] & (~(1<<_bit_index))) | x << _bit_index;
    }
    return *this;
}
Vector<bool>::Reference& Vector<bool>::Reference::operator=( const Vector<bool>::Reference& x ){
    if (_front!=0)
    {
        bool b = (x._front[_int_index] >> x._bit_index) & 1;
        _front[_int_index] = (_front[_int_index] & (~(1<<_bit_index))) | b << _bit_index;
    }
    return *this;
}
Vector<bool>::Reference::operator bool () const{
    if (_front!=0)
    {
        return ((_front[_int_index] >> _bit_index ) & 1);
    } else { 
        return false;
    } 
}
template <typename T>
Vector<T>::Vector() : _size(0), _capacity(2) 
{
    try{
        _front = new T[_capacity]();
    } catch (const std::bad_alloc e) {
        _front = 0;
        _capacity = 0;
        throw;
    }
}
//specialization:
// template<> not used for a member of a specialization
Vector<bool>::Vector() : _size(0) {
    _capacity = _chunk;
    //en 32 lång vektor med 0:or
    _front = new unsigned int[1]();
} 

template <typename T>
Vector<T>::Vector(std::size_t size) : _size(size) 
{
	//capacity är 2 som minst
    _capacity = (max(2, size*1.5));
    //kolla efter int overflow:
    _capacity = (_capacity<=size ? size : _capacity);
    try {
        _front = new T[_capacity]();
    } catch (const std::exception& e) { //generellt exception, både new och T() kan kasta.
        _size = 0;
        _capacity = 0;
        _front = 0;
        throw;
    }
}

Vector<bool>::Vector(std::size_t size){
    std::size_t intern_size = size_convert(size);
    _front = new unsigned int[intern_size]();
    _size = size;
    _capacity = intern_size*_chunk;
}

template <typename T>
Vector<T>::Vector(std::size_t size, T value) : Vector(size)
{
    try{
        for (std::size_t i = 0; i < _size; ++i)
        {
            using std::swap;
            T copy(value); //kan kasta!
            swap(_front[i],copy);
        }
    } catch (const std::exception e){   //om T copy(value) kastar
        //ta bort det som allokerats:
        delete [] _front;
        _front = 0; //så det inte blir dubbel delete av destruktorn.
        _size = 0;
        _capacity = 0;
        throw;
    }
}

Vector<bool>::Vector(std::size_t size, bool value) : Vector(size)
{
    //om value är 0 behöver inget göras. 
    //Men om value är 1 måste intarnas alla bitar sättas till 1.
    std::size_t intern_size = size_convert(size);
    if (value)
    {
        for (std::size_t i = 0; i < intern_size; ++i)
        {
            _front[i] = std::numeric_limits<unsigned int>::max(); //sätt till maximalt värde för unsigned int
        }
    }
}

template <typename T>
Vector<T>::Vector(std::initializer_list<T> il) 
: Vector(il.size())
{
    //lägg över värdena från initializer list:
    try{
        auto iter = _front;
        using std::swap;
        for (T x : il) //kan kasta, copy
        {
            swap(*iter, x);
            ++iter;
        }
    } catch (const std::exception e){
        delete [] _front;
        _front = 0;
        _size = 0;
        _capacity = 0;
        throw;
    }   
}

Vector<bool>::Vector(std::initializer_list<bool> il)
: Vector(il.size())
{
    for(std::size_t i = 0; i< il.size(); ++i){
        //vektorns []-operator används nu:
        (*this)[i] = il.begin()[i];
    }
}

template <typename T>
std::size_t Vector<T>::max(const std::size_t& a, const std::size_t& b) const{
  return (a<b)?b:a;
}

std::size_t Vector<bool>::max(const std::size_t& a, const std::size_t& b) const{
  return (a<b)?b:a;
}

template <typename T>
Vector<T>::Vector(const Vector<T>& rhs)
    : _size(rhs._size),  _capacity(rhs._capacity)
{
    try{
        _front = new T[_capacity]();
    } catch (const std::bad_alloc& e){
        _front = 0;
        _size = 0;
        _capacity = 0;
        throw;
    }

    try{
        for (std::size_t i = 0; i < _size; ++i)
        {
            _front[i] = rhs._front[i]; //kopy assignment kan kasta.
        }
    } catch (const std::exception& e){ // om T:s copy assignment har kastat
        delete [] _front;
        _front = 0;
        _size = 0;
        _capacity = 0;
        throw;
    }   
}

Vector<bool>::Vector(const Vector<bool>& rhs){
    _capacity = rhs._capacity;
    _size = rhs._size;
    std::size_t intern_size = size_convert(_size);
    _front = new unsigned int[intern_size]();
    for (std::size_t i = 0; i < intern_size; ++i)
    {
        _front[i] = rhs._front[i];
    }
}

template <typename T>
Vector<T>& Vector<T>::operator=(const Vector<T>& rhs)
{
    try{
        Vector<T> rhs_copy(rhs); //can throw
        swap(*this, rhs_copy);
    } catch (const std::exception e){
        throw;
    }
    return *this;
}

Vector<bool>& Vector<bool>::operator=(const Vector<bool>& rhs){
    Vector<bool> rhs_copy(rhs);
    swap(*this, rhs_copy);
    return *this;
}

template <typename T>
Vector<T>::Vector(Vector<T>&& rhs)
    : Vector()
{
    swap(*this, rhs); 
    if(this!=&rhs){
        delete [] rhs._front;
        rhs._front = 0;
        rhs._capacity = 0;
    }
}
//TODO:
Vector<bool>::Vector(Vector<bool>&& rhs)
    : Vector()
{
    swap(*this, rhs);
    if(this!=&rhs){
        delete [] rhs._front;
        rhs._front = 0;
        rhs._capacity = 0;
    }
}

template <typename T>
Vector<T>& Vector<T>::operator=(Vector<T>&& rhs)
{
    swap(*this, rhs);
    if(this!=&rhs){
        delete [] rhs._front;
        rhs._front = 0;
        rhs._size = 0;
        rhs._capacity = 0;
    }
    return *this;
}

Vector<bool>& Vector<bool>::operator=(Vector<bool>&& rhs)
{
    swap(*this, rhs);
    if(this!=&rhs){
        delete [] rhs._front;
        rhs._front = 0;
        rhs._size = 0;
        rhs._capacity = 0;
    }
    return *this;
}

template <typename T>
Vector<T>::~Vector(){
    delete [] _front;
    _front = 0;
    _size = 0;
    _capacity = 0;
}
//specialisering
Vector<bool>::~Vector(){
    delete [] _front;
    _front = 0;
    _size = 0;
    _capacity = 0;
}

template <typename T>
void Vector<T>::reallocate()
{
    if(_capacity == std::numeric_limits<std::size_t>::max()) return;

    std::size_t new_cap = max(_capacity, 2)*1.5;
    std::size_t size_copy = _size;
    //check for int overflow:
    new_cap = ( new_cap<_capacity ? std::numeric_limits<std::size_t>::max() : new_cap );
    try{
        Vector<T> copy(new_cap); //can throw
        //kopiera över innehållet till kopian:
        for (std::size_t i = 0; i < _size; ++i)
        {
            T value_copy(_front[i]);
            std::swap(copy[i], value_copy);    //nothrow
        }

        swap(*this, copy);
        _capacity = new_cap;
        _size = size_copy;  //viktigt! blev samma size som cap annars.

    } catch (const std::exception e){
        throw;
    }
}
void Vector<bool>::reallocate(){
    if(_capacity == std::numeric_limits<std::size_t>::max()) return;
    //utöka 
    std::size_t size_copy = _size;
    std::size_t factor = size_convert(_capacity)*1.5 + 0.5;
    std::size_t new_cap = max(_capacity*factor, _chunk);
    //check for int overflow:
    new_cap = ( new_cap<_capacity ? std::numeric_limits<std::size_t>::max() : new_cap );
    Vector<bool> copy(new_cap);

    for (std::size_t i = 0; i < _size; ++i)
        {
            copy[i] = (*this)[i];
        }
    swap(*this, copy);
    //_capacity = new_cap;
    _size = size_copy;
}

template <typename T>
void swap(Vector<T>& first, Vector<T>& second){
    using std::swap; //swap är nothrow
    swap(first._size, second._size);
    swap(first._capacity, second._capacity);
    swap(first._front, second._front);
}

void swap(Vector<bool>& first, Vector<bool>& second){
    using std::swap; //swap är nothrow
    swap(first._size, second._size);
    swap(first._capacity, second._capacity);
    swap(first._front, second._front);
}

template <typename T>
T& Vector<T>::operator[](std::size_t index)
{
    if(index>=_size)
        throw std::out_of_range ("index out of bounds");
    else
        return _front[index];
}
Vector<bool>::Reference Vector<bool>::operator[](std::size_t index)
{
    if(index>=_size)
        throw std::out_of_range ("index out of bounds");
    else{
        Vector<bool>::Reference r;
        r._int_index = intern_indexes(index).first;
        r._bit_index = intern_indexes(index).second;
        r._front = this->_front;
        return r;
    }
}
template <typename T>
const T& Vector<T>::operator[](std::size_t index) const
{
    if(index>=_size)
        throw std::out_of_range ("index out of bounds");
    else
        return _front[index];
}

const Vector<bool>::Reference Vector<bool>::operator[](std::size_t index) const
{
    if(index>=_size)
        throw std::out_of_range ("index out of bounds");
    else{
        Vector<bool>::Reference r;
        r._int_index = intern_indexes(index).first;
        r._bit_index = intern_indexes(index).second;
        r._front = this->_front;
        return r;
    }
}
template <typename T>
void Vector<T>::reset()
{
    try{
        //default -initialisera kopia med copykontruktor:
        Vector<T> default_copy(*this); //can throw
        //reseta alla element:
        for(std::size_t i = 0; i<_size; ++i){
            default_copy[i] = T();
        }
        //byt ut:
        swap(*this, default_copy);
    } catch (const std::exception& e){
        throw;
    }
}
void Vector<bool>::reset()
{
    std::size_t intern_vec_size = size_convert(_capacity);
    for (std::size_t i = 0; i < intern_vec_size; ++i)
    {
        _front[i]=0;
    }
}

//ska tydligen behålla kapaciteten, men storleken ska vara 0.
template <typename T>
void Vector<T>::clear() {
    try{
        Vector<T> copy(*this);
        std::swap(*this, copy);
        _size = 0;
    } catch(const std::exception& e){
        throw;
    }
} 
void Vector<bool>::clear() {
    _size = 0;
}

template <typename T>
std::size_t Vector<T>::size() const {
    return _size;
}

std::size_t Vector<bool>::size() const {
    return _size;
}

template <typename T>
std::size_t Vector<T>::capacity() const {
    return _capacity;
}

std::size_t Vector<bool>::capacity() const {
    return _capacity;
}

template <typename T>
void Vector<T>::push_back(T elem) {
    if (_size>=_capacity)
    {
        reallocate();
    }
    std::swap(_front[_size], elem);
    ++_size;
}

void Vector<bool>::push_back(bool elem) {
    if (_size>=_capacity)
    {
        reallocate();
    }
    //plussar storleken först så att boolvektorns []-operator kan användas utan att kasta exception 
    ++_size;
    (*this)[_size-1] = elem;
}

template <typename T>
void Vector<T>::insert(std::size_t index, T elem) {
    if(index>_size){
        throw std::out_of_range ("index out of bounds");
    }
    else if(index == _size)
    {
        push_back(elem);
    }
    else
    {   
        //first check if vector has enough capacity:
        if(_size>=_capacity)
            reallocate();
        //move elements right
        try{
           //gör kopia av vektorn först:
           Vector<T> tmp_copy(*this);
           //addera ett tomt element i slutet av vektorn, det behövs 
           // annars ger tmp_copy[i],tmp_copy[i-1] out of bound exception
           tmp_copy.push_back(T());
            for (std::size_t i = _size; i > index; --i)
            {
                std::swap(tmp_copy[i],tmp_copy[i-1]);
            }

            std::swap(tmp_copy[index], elem);
            swap(*this, tmp_copy);
            
        } catch (const std::exception& e){
            throw;
        }
    }
}
template <typename T>
void Vector<T>::erase(std::size_t index) {
    if(index>=_size){
        throw std::out_of_range ("index out of bounds");
    } else
    {   
        try{
        //kopia av vektorn:
        Vector<T> tmp_copy(*this); //can throw
        
        //move elements left
        for (std::size_t i = index+1; i < _size; ++i)
        {
            std::swap(tmp_copy[i-1], tmp_copy[i]);
        }
        //byt ut 
        swap(*this, tmp_copy);
    } catch (const std::exception& e){
        throw;
    }
        --_size;
    }
}

template <typename T>
T* Vector<T>::begin() {
    return _front;
}
Vector<bool>::iterator Vector<bool>::begin() {
    iterator mi(*this, 0);
    return mi;
} 

template <typename T>
const T* Vector<T>::begin() const{
    return _front;
}

Vector<bool>::const_iterator Vector<bool>::begin() const {
    const_iterator mci(*this, 0);
    return mci;
}

template <typename T>
T* Vector<T>::end() {
    return _front+_size;
}

Vector<bool>::iterator Vector<bool>::end() {
    iterator mi(*this, _size);
    return mi;
}

template <typename T>
const T* Vector<T>::end() const {
    return _front+_size;
}

Vector<bool>::const_iterator Vector<bool>::end() const {
    const_iterator mci(*this, _size);
    return mci;
}
template <typename T>
T* Vector<T>::find(T const& elem)
{
    T* iter;
    for(iter = _front; iter!=(_front+_size); ++iter )
    {
        if(*iter==elem){
            break; 
        }
    }
    return iter;
}

template <typename T>
const T* Vector<T>::find(T const& elem) const
{
    T* iter;
    for(iter = _front; iter!=(_front+_size); ++iter )
    {
        if(*iter==elem){
            break; 
        }
    }
    return iter;
}

// find iterator
//const iterator:
Vector<bool>::const_iterator::const_iterator() : _vec(0), _index(0){}

Vector<bool>::const_iterator::const_iterator(const Vector<bool>& v, const std::size_t& i) 
: _vec(&v), _index(i) {}

Vector<bool>::const_iterator::const_iterator(const const_iterator& rhs) 
: _vec(rhs._vec), _index(rhs._index) {}

Vector<bool>::const_iterator& Vector<bool>::const_iterator::operator=(const const_iterator& rhs)
{
    _vec = rhs._vec;
    _index = rhs._index;
    return *this;
}

bool operator==(const Vector<bool>::const_iterator& lhs, const Vector<bool>::const_iterator& rhs){
    //måste peka på samma vektor och stå på samma ställe:
    return (lhs._vec == rhs._vec) && (lhs._index == rhs._index); //måste vara samma iterator. Har testat med std::vectors iterator
}
bool operator!=(const Vector<bool>::const_iterator& lhs, const Vector<bool>::const_iterator& rhs){
    return !(lhs==rhs);
}
bool operator<(const Vector<bool>::const_iterator& lhs, const Vector<bool>::const_iterator& rhs){
    return (lhs._index < rhs._index);
}
bool operator<=(const Vector<bool>::const_iterator& lhs, const Vector<bool>::const_iterator& rhs){
    return (lhs._index <= rhs._index);
}
bool operator>(const Vector<bool>::const_iterator& lhs, const Vector<bool>::const_iterator& rhs){
    return (lhs._index > rhs._index);
}
bool operator>=(const Vector<bool>::const_iterator& lhs, const Vector<bool>::const_iterator& rhs){
    return (lhs._index >= rhs._index);
}
std::ptrdiff_t operator-(const Vector<bool>::const_iterator& lhs, const Vector<bool>::const_iterator& rhs){
    return (lhs._index - rhs._index);
}
//std::ptrdiff_t operator-(const Vector<bool>::const_iterator& rhs) const {
//    return _index-rhs._index;
//}

Vector<bool>::const_iterator& Vector<bool>::const_iterator::operator++(){
    ++_index;
    return *this;
}

Vector<bool>::const_iterator& Vector<bool>::const_iterator::operator--(){
    --_index;
    return *this;
}

Vector<bool>::const_iterator Vector<bool>::const_iterator::operator++(int){
    const_iterator ret = *this;
    this->Vector<bool>::const_iterator::operator++();
    return ret;
}

Vector<bool>::const_iterator Vector<bool>::const_iterator::operator--(int){
    const_iterator ret = *this;
    this->Vector<bool>::const_iterator::operator--();
    return ret;
}
Vector<bool>::const_iterator& Vector<bool>::const_iterator::operator+=(int i){
    _index+=i;
    return *this;
}
Vector<bool>::const_iterator& Vector<bool>::const_iterator::operator-=(int i){
    _index-=i;
    return *this;
}

const Vector<bool>::Reference Vector<bool>::const_iterator::operator*() const{
    Vector<bool>::Reference r;
    try{
        r = (*_vec)[_index];
    } catch(const std::out_of_range& e){}
    return r;
}

Vector<bool>::const_iterator Vector<bool>::const_iterator::operator+(int i){
    const_iterator ret = *this;
    ret._index+=i;
    return ret;
}

Vector<bool>::const_iterator Vector<bool>::const_iterator::operator-(int i){
    const_iterator ret = *this;
    ret._index-=i;
    return ret;
}

//non const iterator: ##########################################################
Vector<bool>::iterator::iterator() : _vec(0) {}
Vector<bool>::iterator::iterator(Vector<bool>& v, const std::size_t& i) 
: Vector<bool>::const_iterator(v, i), _vec(&v) {}

Vector<bool>::iterator::iterator(const iterator& rhs) 
: Vector<bool>::const_iterator(rhs), _vec(rhs._vec) {}

Vector<bool>::iterator& Vector<bool>::iterator::operator=(const iterator& rhs)
{
    _vec = rhs._vec;
    _index = rhs._index;
    return *this;
}
Vector<bool>::iterator& Vector<bool>::iterator::operator++(){
    ++_index;
    return *this;
}

Vector<bool>::iterator Vector<bool>::iterator::operator++(int){
    iterator ret = *this;
    this->Vector<bool>::const_iterator::operator++();
    return ret;
}

Vector<bool>::iterator& Vector<bool>::iterator::operator--(){
    --_index;
    return *this;
}

Vector<bool>::iterator Vector<bool>::iterator::operator--(int){
    iterator ret = *this;
    this->Vector<bool>::const_iterator::operator--();
    return ret;
}
Vector<bool>::iterator Vector<bool>::iterator::operator+(int i){
    iterator ret = *this;
    ret._index+=i;
    return ret;
}

Vector<bool>::iterator Vector<bool>::iterator::operator-(int i){
    iterator ret = *this;
    ret._index-=i;
    return ret;
}

Vector<bool>::Reference Vector<bool>::iterator::operator*(){
    try{
        Vector<bool>::Reference r = (*_vec)[_index];
        return r;
    } catch(const std::out_of_range& e){
        Reference r;//dummyreferens, pekarn'te inte på nöt'
        return r;
    }
}
}