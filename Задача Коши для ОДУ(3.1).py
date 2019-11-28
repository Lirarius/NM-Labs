#Задача Коши для обыкновенных дифференциальных уравнений

#№3.1

import numpy as np
import math
import matplotlib.pyplot as plt

#решаем задачу коши y' = f(x, y), y(x0) = y0, x0 <= x <= T для системы методом Рунге—Кутта c шагом h
def runge_kutta(f, x0, y0, T, h):
    xx = [x0]
    yy = [y0]
    
    while(True):
        #предыдущие приближения
        yn = yy[-1]
        xn = xx[-1]
        #заканчиваем, если выходим за пределы отрезка
        if (xn + h > T): 
            break
        #проводим расчет
        k1 = f(xn, yn)
        k2 = f(xn + 0.5 * h, yn + 0.5 * h * k1)
        k3 = f(xn + 0.5 * h, yn + 0.5 * h * k2)
        k4 = f(xn + h, yn + h * k3)
        yn1 = yn + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        #добавляем в массив новые значения
        xx.append(xn + h)
        yy.append(yn1)
    
    return xx, yy

#имеем уравнение:
# u'' = -sin u
#преобразуем к системе:
# u' = v
# v' = -sin u
# u(0) = 1
# v(0) = u'(0) = 0

#решаем численно
xx, uu = runge_kutta(lambda x, uv: np.array([uv[1], -math.sin(uv[0])]), 0, np.array([1, 0]), 4 * math.pi, 0.01)
#uu - содержит вектора, которые содержат функцию и ее производную
#выберем только функцию
uu = [x[0] for x in uu]

#построим график
plt.plot(xx, uu)
#сохраним график в файл runge_demo_result.png
plt.savefig('runge_demo_result.png')