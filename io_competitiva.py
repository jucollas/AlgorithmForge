# ===============================================
# ENTRADA Y SALIDA EN PYTHON (PROGRAMACION COMPETITIVA)
# ===============================================
# Autor: Ejemplo generado por ChatGPT
# Descripcion:
# Este archivo muestra distintas maneras eficientes de
# recibir e imprimir datos en Python para programacion competitiva.

import sys

# ------------------------------
# 1. Lectura basica de una linea
# ------------------------------
# sys.stdin.readline() lee una linea completa mas rapido que input().
linea = sys.stdin.readline().strip()
print("Entrada leida:", linea)

# ------------------------------
# 2. Leer varios numeros en una linea
# ------------------------------
a, b = map(int, sys.stdin.readline().split())
print("Suma:", a + b)

# ------------------------------
# 3. Leer multiples lineas con cantidad conocida
# ------------------------------
n = int(sys.stdin.readline())
for i in range(n):
    x, y = map(int, sys.stdin.readline().split())
    print(f"Resultado {i+1}:", x + y)

# ------------------------------
# 4. Redefinir input para simplicidad
# ------------------------------
input = sys.stdin.readline  # ahora input() sera rapido
n = int(input())
for _ in range(n):
    a, b = map(int, input().split())
    print(a * b)

# ------------------------------
# 5. Salida rapida con sys.stdout.write
# ------------------------------
sys.stdout.write("Esto es una salida rapida\n")

# ------------------------------
# 6. Salida eficiente acumulando resultados
# ------------------------------
input = sys.stdin.readline
n = int(input())
resultados = []
for _ in range(n):
    a, b = map(int, input().split())
    resultados.append(str(a - b))
sys.stdout.write("\n".join(resultados) + "\n")

# ------------------------------
# 7. Leer hasta EOF (fin de archivo)
# ------------------------------
# Muy util cuando no se da el numero de casos.
for linea in sys.stdin:
    datos = linea.strip().split()
    if not datos:
        continue
    nums = list(map(int, datos))
    print("Suma de la linea:", sum(nums))
