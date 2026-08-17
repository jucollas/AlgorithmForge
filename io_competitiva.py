import sys
sys.set_int_max_str_digits(10**9)
sys.setrecursionlimit(10**9)

input = sys.stdin.readline  # ahora input() sera rapido

# ------------------------------
# 1. Lectura basica de una linea
# ------------------------------
# sys.stdin.readline() lee una linea completa mas rapido que input().
linea = sys.stdin.readline().strip()

# ------------------------------
# 2. Leer varios numeros en una linea
# ------------------------------
a, b = map(int, sys.stdin.readline().split())

# ------------------------------
# 3. Leer multiples lineas con cantidad conocida
# ------------------------------
n = int(sys.stdin.readline())
for i in range(n):
    x, y = map(int, sys.stdin.readline().split())

# ------------------------------
# 5. Salida rapida con sys.stdout.write
# ------------------------------
sys.stdout.write("Esto es una salida rapida\n")

# ------------------------------
# 7. Leer hasta EOF (fin de archivo)
# ------------------------------
# Muy util cuando no se da el numero de casos.
for linea in sys.stdin:
    datos = linea.strip().split()
    if not datos:
        continue
    nums = list(map(int, datos))
