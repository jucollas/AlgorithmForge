### returns the shift giving
### lexicographically minimal string
### among trans. phi(c_1S) |-> Sc_1
def booth(s: str) -> int:
    s = s + s
    n = len(s) // 2
    i, j, k = 0, 1, 0

    while i < n and j < n and k < n:
        if s[i + k] == s[j + k]:
            k += 1
        elif s[i + k] > s[j + k]:
            i = i + k + 1
            if i <= j:
                i = j + 1
            k = 0
        else:
            j = j + k + 1
            if j <= i:
                j = i + 1
            k = 0
    return min(i, j)