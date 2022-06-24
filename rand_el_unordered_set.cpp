/*
1. check the unordered set in not empty (if it is, there is no hope)
2. generate a random value_type element
3. if already in the unordered set return it else insert it
4. get an iterator it on this element
5. get the random element as *(it++) (and if *it is the last element the get the first element)
6. delete the element you inserted and return the value in (5)
*/

#include <iostream>
#include <unordered_set>
using namespace std;

int random_element_set(unordered_set<int> &s, int r) {
    if(s.empty()){
        printf("PRAZAN SKUP!");
        return -1;
    }
    if (s.find(r) != s.end()) {
        return r;
    } else {
        s.insert(r);
    }
    auto it = s.find(r);
    advance(it, 1);
    if (it == s.end()) {
        s.erase(r);
        return *(s.begin());
    } else {
        s.erase(r);
        return *(it);
    }
}