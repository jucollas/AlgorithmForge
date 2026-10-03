/*
Author: Juan Diego Collazos
Date: Fri Oct  2 19:30:31 2026
*/

/*
The Knuth-Morris-Pratt algorithm
Search for a pattern inside a text.
pattern = string we are searching for
text    = string where we search
O(|pattern| + |text|)
*/
vector<int> kmp(const string& pattern, const string& text) {
    vector<int> ans;
    if (pattern.empty()) return ans;
    vector<int> pi = prefix_function(pattern);
    int j = 0;
    for (int i = 0; i < sz(text); ++i) {
        while (j > 0 && text[i] != pattern[j]) j = pi[j - 1];
        if (text[i] == pattern[j]) ++j;
        if (j == sz(pattern)) {
            ans.push_back(i + 1 - sz(pattern));
            j = pi[j - 1];
        }
    }
    return ans;
}

/*
Counting the number of occurrences of each prefix of s
inside t.
ans[i] = number of occurrences of
         s[0 ... i-1] inside t.
O(|s| + |t|)
*/
vector<int> count_occ_each_prefix(
    const string& s,
    const string& t
) {
    const int n = sz(s);
    vector<int> ans(n + 1);
    if (n == 0) return ans;
    vector<int> pi = prefix_function(s);
    int j = 0;
    for (int i = 0; i < sz(t); ++i) {
        while (j > 0 && t[i] != s[j]) j = pi[j - 1];
        if (t[i] == s[j]) ++j;
        ++ans[j];
        if (j == n) j = pi[j - 1];
    }
    for (int i = n; i > 0; --i) ans[pi[i - 1]] += ans[i];
    return ans;
}

/*
Counting the number of occurrences of each prefix
inside the string itself.
ans[i] = number of occurrences of
         s[0 ... i-1] inside s.
O(n)
*/
vector<int> count_occ_each_prefix_itself(const string& s) {
    const int n = sz(s);
    vector<int> ans(n + 1);
    if (n == 0) return ans;
    vector<int> pi = prefix_function(s);
    /*Each pi[i] tells us that the prefix of length pi[i]
    occurs as a suffix ending at position i.*/
    for (int i = 0; i < n; ++i) ++ans[pi[i]];
    /* Propagate occurrences through the border tree. */
    for (int i = n; i > 0; --i) ans[pi[i - 1]] += ans[i];
    /* Every prefix occurs at least once: as itself.*/
    for (int i = 0; i <= n; ++i) ++ans[i];

    return ans;
}

/*
The number of different substrings in a string.
For every prefix s[0...i]:
    new substrings = length - maximum prefix-function value
                     of reverse(s[0...i])
O(n²)
*/
lint number_diff_substring(const string& s) {
    const int n = sz(s);
    lint ans = 0;
    string t;
    vector<int> pi;
    for (int i = 0; i < n; ++i) {
        t += s[i];
        string tinv(t.rbegin(), t.rend());
        pi = prefix_function(tinv);
        int mx = *max_element(all(pi));
        ans += sz(t) - mx;
    }
    return ans;
}

/*
Compressing a string.
Returns the shortest string whose repetition
produces s.
O(n)
*/
string compress_string(const string& s) {
    const int n = sz(s);
    if (n == 0) return "";
    vector<int> pi = prefix_function(s);
    int k = n - pi[n - 1];
    if (n % k == 0) return s.substr(0, k);
    return s;
}
