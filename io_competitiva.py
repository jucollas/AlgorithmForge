## necesities
import sys
sys.set_int_max_str_digits(10**9)
sys.setrecursionlimit(10**9)
## Examples
linea = sys.stdin.readline().strip()
arr = list(map(int, sys.stdin.readline().split()))
sys.stdout.write("Esto es una salida rapida\n")
## leer hasta EOF
for linea in sys.stdin:
    datos = linea.strip().split()
    if not datos: continue
    nums = list(map(int, datos))