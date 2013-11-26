#ifndef RANGE_HPP
#define RANGE_HPP

//This is a light version of boost::interval
//Most code is copied from boost::interval library

#include <algorithm>

template<class T>
class range_base
{
public:
  T const &lower() const;
  T const &upper() const;

  range_base();
  range_base(T const &v);
  range_base(T const &l, T const &u);
  range_base(range_base<T> const &r);

  range_base& operator=(range_base<T> const &r);


private:
  T low;
  T up;
  
};

template<class T> inline
bool empty(const range_base<T>& x)
{
  return x.upper() < x.lower();
}

template<class T> inline
range_base<T> intersect(const range_base<T>& x,
			const range_base<T>& y)
{
  return range_base<T>(std::max(x.lower(),y.lower()),
		       std::min(x.upper(),y.upper()));
}

template<class T> inline
range_base<T> operator-(const range_base<T>& x,
			const T&y)
{
  return range_base<T>(x.lower()-y, x.upper()-y);
}

template<class T> inline
range_base<T> operator-(const T& x,
			const range_base<T>&y)
{
  return range_base<T>(x-y.upper(), x-y.lower());
}


//y MUST be non-negtive
template<class T> inline
range_base<T> operator*(const range_base<T>& x,
			const T&y)
{
  return range_base<T>(x.lower()*y, x.upper()*y);
}
//x MUST be non-negtive
template<class T> inline
range_base<T> operator*(const T& x,
			const range_base<T>&y)
{
  return range_base<T>(x*y.lower(), x*y.upper());
}



template<class T> inline
const T& range_base<T>::upper() const
{
  return up;
}

template<class T> inline
T upper(range_base<T> const &v)
{
  return v.upper();
}


template<class T> inline
const T& range_base<T>::lower() const
{
  return low;
}

template<class T> inline
T lower(range_base<T> const &v)
{
  return v.lower();
}



template<class T> inline
range_base<T>::range_base():
  low(static_cast<T>(0)), up(static_cast<T>(0))
{}

template<class T> inline
range_base<T>::range_base(T const &v): low(v), up(v)
{}

template<class T> inline
range_base<T>::range_base(T const &l, T const &h): low(l), up(h)
{}


template<class T> inline
range_base<T>::range_base(range_base<T> const &r):low(r.lower()),up(r.upper())
{}


template<class T> inline
range_base<T>& range_base<T>::operator=(range_base<T> const &r)
{
  low = r.lower();
  up  = r.upper();
  return *this;
}


#endif
