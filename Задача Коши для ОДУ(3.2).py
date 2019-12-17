#Метод предиктор–корректор

import numpy as np

def predict_correct(f, x0, y0, T, h):

    xx = [x0]
    yy = [y0]
    
    while True:
        if xx[-1] + h > T: 
            break
            
        if len(yy) < 4:
            yy.append(yy[-1] + h * f(xx[-1], yy[-1]))
            xx.append(xx[-1] + h)
        else:
            yn = yy[-1] + h * (55 * f(xx[-1], yy[-1]) - 58 * f(xx[-2], yy[-2]) + 37 * f(xx[-3], yy[-3]) - 9 * f(xx[-4], yy[-4])) / 24.0
            yn = yy[-1] + h * (9 * f(xx[-1] + h, yn) + 19 * f(xx[-1], yy[-1]) - 5 * f(xx[-2], yy[-2]) + f(xx[-3], yy[-3])) / 24.0
            
            yy.append(yn)
            xx.append(xx[-1] + h)
    return xx, yy
    
a = 0.2
b = 0.2
c = 2.5
xx, uu = predict_correct(lambda x, uv: np.array([uv[1] - uv[2], uv[0] + a * uv[1], b + uv[2] * (uv[0] - c)]), 0, np.array([1, 1, 1]), 10, 0.010)

for i in range(len(xx)):
    print("u1(" + str(xx[i]) + ") = " + str(uu[i][0]))
for i in range(len(xx)):
    print("u2(" + str(xx[i]) + ") = " + str(uu[i][1]))
for i in range(len(xx)):
    print("u3(" + str(xx[i]) + ") = " + str(uu[i][2]))

a = 0.2
b = 0.2
c = 5.0
xx, uu = predict_correct(lambda x, uv: np.array([uv[1] - uv[2], uv[0] + a * uv[1], b + uv[2] * (uv[0] - c)]), 0, np.array([1, 1, 1]), 10, 0.010)

for i in range(len(xx)):
    print("u1(" + str(xx[i]) + ") = " + str(uu[i][0]))
for i in range(len(xx)):
    print("u2(" + str(xx[i]) + ") = " + str(uu[i][1]))
for i in range(len(xx)):
    print("u3(" + str(xx[i]) + ") = " + str(uu[i][2]))