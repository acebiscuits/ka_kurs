# -*- coding: utf-8 -*-
from sage.all import *
import numpy as np
from math import log, floor
from random import randint

def LpCalculating(p):
    log_p = math.log(p)
    log_log_p = math.log(log_p)
    return math.exp((1 / math.sqrt(2)) * (log_p**0.5) * (log_log_p**0.5))

def discrCalculating(a, b, N):
    return (4*a**3 + 27*b**2) % N

def binaryDecomposition(k):
    bits = bin(k)[2:]

    positions = [i for i, bit in enumerate(reversed(bits)) if bit == '1']

    return bits, positions

def inverseMod(n, N):
    try:
        return pow(n, -1, N)
    except ValueError:
        return None
    
def pointDouble(x1, y1, A, N):
    if y1 == 0:
        return None

    numerator = (3 * x1**2 + A) % N
    denominator = (2 * y1) % N

    invDenominator = inverseMod(denominator, N)
    if invDenominator is None:
        d = math.gcd(denominator, N)
        if 1 < d < N:
            raise ZeroDivisionError(f"Делитель найден: {d}")
        else:
            raise ZeroDivisionError("Обратного элемента не нашлось и делитель не найден, надо пробовать другую уривую")

    lambd = (numerator * invDenominator) % N

    x3 = (lambd**2 - 2 * x1) % N
    y3 = (lambd * (x1 - x3) - y1) % N

    return x3, y3

def pointAddition(x1, y1, x2, y2, A, N):
    if x1 is None and y1 is None:#P+Q=Q
        return x2, y2
    if x2 is None and y2 is None:#P+Q=P
        return x1, y1

    if x1 == x2 and (y1 + y2) % N == 0:#P+(−P)=O
        return None, None

    if x1 == x2 and y1 == y2:
        # Удвоение точки
        return pointDouble(x1, y1, A, N)
    else:
        # Сложение разных точек
        numerator = (y2 - y1) % N
        denominator = (x2 - x1) % N

    invDenominator = inverseMod(denominator, N)
    if invDenominator is None:
        d = math.gcd(denominator, N)
        if 1 < d < N:
            raise ZeroDivisionError(f"Делитель найден: {d}")
        else:
            raise ZeroDivisionError("Обратного элемента нет и делитель не найден")

    lambd = (numerator * invDenominator) % N

    x3 = (lambd**2 - x1 - x2) % N
    y3 = (lambd * (x1 - x3) - y1) % N

    return x3, y3

def scalarMultiplication(k, x, y, A, N):
    R = (None, None)  # нейтральный элемент (точка на бесконечности)
    Q = (x, y)        # исходная точка P

    bits = bin(k)[2:]  # двоичное представление k, строка без '0b'

    for bit in bits:
        Q = pointDouble(Q[0], Q[1], A, N)

        if bit == '1':
            R = pointAddition(R[0], R[1], Q[0], Q[1], A, N)

    return R


def lenstraEcm(N, B1, B2, tries):
    for tryNumber in range(1, tries):
        x = randint(0, N - 1)
        y = randint(0, N - 1)
        a = randint(0, N - 1)

        b = (y^2 - x^3 - a*x) % N

        discr = discrCalculating(a, b, N)
        
        gcd = math.gcd(discr, N)
        if gcd == N:
            continue
        if 1 < gcd < N:
            return gcd
        
        for p in prime_range(2, B1):
            e = floor(log(B2) / log(p))
            if e == 0:
                continue
            k = p**e

            try:
                binaryDecomposition(k)
                
                R = scalarMultiplication(k, x, y, a, N)
                #xres, yres = R[0], R[1]

            except ZeroDivisionError as e:
                msg = str(e)
                if "Делитель найден" in msg:
                    d = int(msg.split(":")[-1].strip())
                    return d
                else:
                    continue

            except (ArithmeticError, ValueError) as err:
                print("Ошибка арифметики:", err)

    print("Не удалось найти делитель.")
    return None


bits_p = 13
bits_q = 30

p = random_prime(2**bits_p, lbound=2**(bits_p - 1))
print(f"p = {p} (битность: {p.nbits()})")
q = random_prime(2**bits_q, lbound=2**(bits_q - 1))
print(f"q = {q} (битность: {q.nbits()})")

N = p * q
print(f"битность N ({N}): {N.nbits()})")
print(f"Ищем делитель числа N = {N} = {p} * {q}")

Lp = LpCalculating(p)
print(f"L_p({p}) = {Lp}")
B1 = int(Lp)
print(f"B1 = {B1}")
mulNumber = 5
B2 = mulNumber * B1
print(f"B2 = {B2}")
tries = 200000

delimeter = lenstraEcm(N, B1, B2, tries)

if delimeter:
    print(f"Делитель найден: {delimeter}, второй: {N // delimeter}")
else:
    print("Делитель не найден")
