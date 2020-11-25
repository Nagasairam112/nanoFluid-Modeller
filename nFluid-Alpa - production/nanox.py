import os
import sys
import pandas as pd
from sklearn.metrics import r2_score
import numpy as np
from PyQt5 import QtWidgets, uic
from PyQt5.QtWidgets import QFileDialog, QMessageBox
qtcreator_file = "nanoflui2.ui"
Ui_MainWindow, QtBaseClass = uic.loadUiType(qtcreator_file)
class MyWindow(QtWidgets.QMainWindow, Ui_MainWindow):
    def __init__(self):
        QtWidgets.QMainWindow.__init__(self)
        Ui_MainWindow.__init__(self)
        self.setupUi(self)
        self.tabWidget.setTabVisible(0, True)
        self.tabWidget.setTabVisible(1, False)
        self.tabWidget.setTabVisible(2, False)
        self.tabWidget.setTabVisible(3, False)
        self.tabWidget.setTabVisible(4, False)
        self.tabWidget.setTabVisible(5, False)

        self.sai.setVisible(False)
        self.jbox.setVisible(False)
        self.z3.toggled.connect(self.showbox)
        self.t1.toggled.connect(self.hx)
        self.t2.toggled.connect(self.hx)
        self.t3.toggled.connect(self.sx)
        self.jo2.toggled.connect(self.sf)
        self.sav.clicked.connect(self.save)
        self.sh.clicked.connect(self.tudf)
        self.bill.setVisible(False)
        self.bill.clicked.connect(self.tdef)

        self.button1.clicked.connect(self.viz)
        self.button2.clicked.connect(self.print)
        self.button3.clicked.connect(self.save)
        self.button5.clicked.connect(self.texs)
        self.button6.clicked.connect(self.ted)
        self.button9.clicked.connect(self.exp)
        self.so.clicked.connect(self.rowadd)
        self.jup.clicked.connect(self.act)
        self.about.clicked.connect(self.abo)
        self.button10.clicked.connect(self.ro)
        self.button11.clicked.connect(self.hyp)
        self.so_2.clicked.connect(self.adr)
        self.button12.clicked.connect(self.tphi)
        self.button3.setEnabled(False)


    def tdef(self):
        try:
            self.MplWidget2.canvas.Axes3D.clear()
            self.MplWidget2.canvas.Axes3D2.clear()
            self.MplWidget2.canvas.Axes3D3.clear()
            self.MplWidget2.canvas.Axes3D4.clear()

            self.MplWidget2.canvas.Axes3D.grid(True)
            self.MplWidget2.canvas.Axes3D2.grid(True)
            self.MplWidget2.canvas.Axes3D3.grid(True)
            self.MplWidget2.canvas.Axes3D4.grid(True)

            self.MplWidget2.canvas.Axes3D.set_xlabel("Temperature (C)")
            self.MplWidget2.canvas.Axes3D.set_ylabel("Concentration (%)")
            self.MplWidget2.canvas.Axes3D.set_zlabel("Density (Kg/m^3)")

            self.MplWidget2.canvas.Axes3D2.set_xlabel("Temperature (C)")
            self.MplWidget2.canvas.Axes3D2.set_ylabel("Concentration (%)")
            self.MplWidget2.canvas.Axes3D2.set_zlabel("Specific heat (Kg/m.K)")

            self.MplWidget2.canvas.Axes3D3.set_xlabel("Temperature (C)")
            self.MplWidget2.canvas.Axes3D3.set_ylabel("Concentration (%)")
            self.MplWidget2.canvas.Axes3D3.set_zlabel("Thermal Conductivity (W/m.K)")

            self.MplWidget2.canvas.Axes3D4.set_xlabel("Temperature (C)")
            self.MplWidget2.canvas.Axes3D4.set_ylabel("Concentration (%)")
            self.MplWidget2.canvas.Axes3D4.set_zlabel("Viscoscity (Kg/m.s)")
            self.MplWidget2.canvas.draw()

            phi = float(self.sika2.text())
            tem = float(self.v8.text())
            co = np.linspace(0.0, phi / 100, num=11)
            temp = np.linspace(0.0, tem, num=11)
            t = (self.nams.text())
            c = str(t)
            self.MplWidget2.canvas.figure.suptitle("Thermo Physical Properties of " + c)
            model = self.table.model()
            model2 = self.table3.model()

            val = model.index(0, 0)
            val1 = model.index(0, 1)
            val2 = model.index(0, 2)
            val3 = model.index(0, 3)

            va = model.index(1, 0)
            va1 = model.index(1, 1)
            va2 = model.index(1, 2)
            va3 = model.index(1, 3)

            v = model.index(2, 0)
            v1 = model.index(2, 1)
            v2 = model.index(2, 2)
            v3 = model.index(2, 3)

            p1 = float(model.data(val))
            p2 = float(model.data(val1))
            p3 = float(model.data(val2))
            p4 = float(model.data(val3))

            c1 = float(model.data(va))
            c2 = float(model.data(va1))
            c3 = float(model.data(va2))
            c4 = float(model.data(va3))

            b1 = float(model.data(v))
            b2 = float(model.data(v1))
            b3 = float(model.data(v2))
            b4 = float(model.data(v3))

            d_n = p1 + (p2) * temp + (p3) * pow(temp, 2) + (p4) * pow(temp, 3)
            c_p = c1 + (c2) * temp + (c3) * pow(temp, 2) + (c4) * pow(temp, 3)
            k_p = b1 + (b2) * temp + (b3) * pow(temp, 2) + (b4) * pow(temp, 3)
            j1 = float(self.x1.text())
            j2 = float(self.x2.text())
            j3 = float(self.x3.text())

            q = model2.index(0, 0)
            q1 = model2.index(0, 1)
            q2 = model2.index(0, 2)
            q3 = model2.index(0, 3)

            a = model2.index(1, 0)
            a1 = model2.index(1, 1)
            a2 = model2.index(1, 2)
            a3 = model2.index(1, 3)

            ve = model2.index(2, 0)
            ve1 = model2.index(2, 1)
            ve2 = model2.index(2, 2)
            ve3 = model2.index(2, 3)

            wo = model2.index(3, 0)
            wo1 = model2.index(3, 1)
            wo2 = model2.index(3, 2)
            wo3 = model2.index(3, 3)

            c1 = float(model2.data(q))
            c2 = float(model2.data(q1))
            c3 = float(model2.data(q2))
            c4 = float(model2.data(q3))

            k1 = float(model2.data(a))
            k2 = float(model2.data(a1))
            k3 = float(model2.data(a2))
            k4 = float(model2.data(a3))

            l1 = float(model2.data(ve))
            l2 = float(model2.data(ve1))
            l3 = float(model2.data(ve2))
            l4 = float(model2.data(ve3))

            m1 = float(model2.data(wo))
            m2 = float(model2.data(wo1))
            m3 = float(model2.data(wo2))
            m4 = float(model2.data(wo3))

            d_w = (c1) + (c2) * temp + (c3) * pow(temp, 2) + (c4) * pow(temp, 3)
            k_w = (k1) + (k2) * temp + (k3) * pow(temp, 2) + (k4) * pow(temp, 3)
            mu_w = (l1) + (l2) * temp + (l3) * pow(temp, 2) + (l4) * pow(temp, 3)
            c_w = (m1) + (m2) * temp + (m3) * pow(temp, 2) + (m4) * pow(temp, 3)

            d_nf = (1 - co) * d_w + co * d_n
            cp_nf = (1 - co) * c_w + co * c_p
            k_nf = ((k_p + 2 * k_w + 2 * co * (k_p - k_w)) / (k_p + 2 * k_w - co * (k_p - k_w))) * k_w
            mu_nf = mu_w * (j1 + j2 * phi + j3 * pow(phi, 2))

            self.MplWidget2.canvas.Axes3D.plot(temp, co, d_nf, 'darkorange')
            self.MplWidget2.canvas.Axes3D2.plot(temp, co, cp_nf, 'darkgreen')
            self.MplWidget2.canvas.Axes3D3.plot(temp, co, k_nf, 'darkblue')
            self.MplWidget2.canvas.Axes3D4.plot(temp, co, mu_nf, 'darkviolet')
            self.MplWidget2.canvas.draw()
        except ValueError:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Warning)
            msg.setWindowTitle("Warning")
            msg.setText("Chek the input format")

            def msgButtonClick(i):
                msg.close()

            msg.setStandardButtons(QMessageBox.Cancel)
            msg.show()
            msg.buttonClicked.connect(msgButtonClick)


    def tudf(self):
        try:
         self.tex.setText("")
         x=' '
         model5 = self.huc1.model()
         model6 = self.huc2.model()
         print("1")
         pho = float(self.aha2.text())
         te = (self.aha1.text())
         con = str(te)
         print(pho,con)
         q1 = model5.index(0, 0)
         q2 = model5.index(0, 1)
         q3 = model5.index(0, 2)
         q4 = model5.index(0, 3)

         a = model5.index(1, 0)
         a1 = model5.index(1, 1)
         a2 = model5.index(1, 2)
         a3 = model5.index(1, 3)

         ve = model5.index(2, 0)
         ve1 = model5.index(2, 1)
         ve2 = model5.index(2, 2)
         ve3 = model5.index(2, 3)

         wo = model5.index(3, 0)
         wo1 = model5.index(3, 1)
         wo2 = model5.index(3, 2)
         wo3 = model5.index(3, 3)

         e1 = float(model5.data(q1))
         e2 = float(model5.data(q2))
         e3 = float(model5.data(q3))
         e4 = float(model5.data(q4))

         k1 = float(model5.data(a))
         k2 = float(model5.data(a1))
         k3 = float(model5.data(a2))
         k4 = float(model5.data(a3))

         l1 = float(model5.data(ve))
         l2 = float(model5.data(ve1))
         l3 = float(model5.data(ve2))
         l4 = float(model5.data(ve3))

         m1 = float(model5.data(wo))
         m2 = float(model5.data(wo1))
         m3 = float(model5.data(wo2))
         m4 = float(model5.data(wo3))

         val = model6.index(0, 0)
         val1 = model6.index(0, 1)
         val2 = model6.index(0, 2)
         val3 = model6.index(0, 3)

         va = model6.index(1, 0)
         va1 = model6.index(1, 1)
         va2 = model6.index(1, 2)
         va3 = model6.index(1, 3)

         v = model6.index(2, 0)
         v1 = model6.index(2, 1)
         v2 = model6.index(2, 2)
         v3 = model6.index(2, 3)

         p1 = float(model6.data(val))
         p2 = float(model6.data(val1))
         p3 = float(model6.data(val2))
         p4 = float(model6.data(val3))

         c1 = float(model6.data(va))
         c2 = float(model6.data(va1))
         c3 = float(model6.data(va2))
         c4 = float(model6.data(va3))

         b1 = float(model6.data(v))
         b2 = float(model6.data(v1))
         b3 = float(model6.data(v2))
         b4 = float(model6.data(v3))

         self.tex.append('/****************ANSYS FLUENT UDF FOR ' + (con) + +1 * x + 'NANOFLUID****************/')
         self.tex.append('/****Compile before starting the Fluent****/')
         self.tex.append('#include "udf.h"')
         self.tex.append('#define phi %f ' % (pho / 100))
         self.tex.append('DEFINE_PROPERTY(nanofluid_conductivity,c,t)')
         self.tex.append('{')
         self.tex.append('real k_w;')
         self.tex.append('real k_nf;')
         self.tex.append('real ctemp = C_T(c,t);')
         self.tex.append('k_w=%f+(%f)*ctemp+(%f)*ctemp*ctemp+(%f)*ctemp*ctemp*ctemp;' % (l1, l2, l3, l4))
         self.tex.append('k_p=%f+%f*ctemp+(%f)*ctemp*ctemp+(%f)*ctemp*ctemp*ctemp;'% (b1, b2, b3, b4))
         self.tex.append('k_nf=(((k_p+(2*k_w)-(2*(k_w-k_p))*phi)/(k_p+(2*k_w)+(k_w-k_p)*phi)))*k_w;')
         self.tex.append('return k_nf;')
         self.tex.append('}')
         self.tex.append('DEFINE_PROPERTY(nanofluid_viscosity,c,t)')
         self.tex.append('{')
         self.tex.append('real mu_w;')
         self.tex.append('real mu_nf;')
         self.tex.append('real ctemp2 = C_T(c,t);')
         self.tex.append('mu_nf = %f+(%f)*ctemp+(%f)*ctemp*ctemp+(%f)*ctemp*ctemp*ctemp' % (m1, m2, m3, m4))
         self.tex.append('return mu_nf;')
         self.tex.append('}')
         self.tex.append('DEFINE_SPECIFIC_HEAT(nanofluid_specificheat,c,t)')
         self.tex.append('{')
         self.tex.append('real cp_w;')
         self.tex.append('real cp_nf;')
         self.tex.append('real ctemp3 = C_T(c,t);')
         self.tex.append('cp_w=%f+(%f)*ctemp+(%f)*ctemp*ctemp+(%f)*ctemp*ctemp*ctemp;' % (k1, k2, k3, k4))
         self.tex.append('cpp=%f+(%f)*ctemp+(%f)*ctemp*ctemp+(%f)*ctemp*ctemp*ctemp;' % (c1, c2, c3, c4))
         self.tex.append('cp_nf = (1-phi)*cp_w+phi*cpp;')
         self.tex.append('return cp_nf;')
         self.tex.append('}')
         self.tex.append('DEFINE_PROPERTY(nanofluid_density,c,t)')
         self.tex.append('{')
         self.tex.append('real rho_w;')
         self.tex.append('real rho_nf;')
         self.tex.append('real ctemp4 = C_T(c,t);')
         self.tex.append('rho_w = %f+(%f)*ctemp+(%f)*ctemp*ctemp+(%f)*ctemp*ctemp*ctemp;' % (e1, e2, e3, e4))
         self.tex.append('rho_p = %f+(%f)*ctemp+(%f)*ctemp*ctemp+(%f)*ctemp*ctemp*ctemp;' % (p1, p2, p3, p4))
         self.tex.append('rho_nf = (1-phi)*rho_w+phi*rho_p;')
         self.tex.append('}')
         self.sav.setVisible(True)


        except ValueError:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Warning)
            msg.setWindowTitle("Warning")
            msg.setText("Input Format must be a Numeric Values")

            def msgButtonClick(i):
                msg.close()

            msg.setStandardButtons(QMessageBox.Cancel)
            msg.show()
            msg.buttonClicked.connect(msgButtonClick)


    def sf(self):
        if self.jo2.isChecked()==True:
         self.jbox.setVisible(True)
         self.tex.setText("")
         self.button3.setVisible(False)
         self.sav.setVisible(False)
        else:
         self.jbox.setVisible(True)
         self.button3.setVisible(True)
         self.print()


    def tphi(self):
        try:
            self.MplWdiget9.canvas.A.clear()
            self.MplWdiget9.canvas.A2.clear()
            self.MplWdiget9.canvas.A3.clear()
            self.MplWdiget9.canvas.A4.clear()

            self.MplWdiget9.canvas.A.set_xlabel('phi/T')
            self.MplWdiget9.canvas.A.set_ylabel('Density')
            self.MplWdiget9.canvas.A2.set_xlabel('phi/T')
            self.MplWdiget9.canvas.A2.set_ylabel('Thermal Conductivity')
            self.MplWdiget9.canvas.A3.set_xlabel('phi/T')
            self.MplWdiget9.canvas.A3.set_ylabel('Specific Heat')
            self.MplWdiget9.canvas.A4.set_xlabel('phi/T')
            self.MplWdiget9.canvas.A4.set_ylabel('Dynamic Viscocity')

            self.MplWdiget9.canvas.A.grid(True)
            self.MplWdiget9.canvas.A2.grid(True)
            self.MplWdiget9.canvas.A3.grid(True)
            self.MplWdiget9.canvas.A4.grid(True)
            self.MplWdiget9.canvas.draw()

            col_count = self.table5.columnCount()
            row_count = self.table5.rowCount()
            headers = [str(self.table5.horizontalHeaderItem(i).text()) for i in range(col_count)]
            df_list = []
            for row in range(row_count):
                df_list2 = []
                for col in range(col_count):
                    table_item = self.table5.item(row, col)
                    df_list2.append('' if table_item is None else str(table_item.text()))
                df_list.append(df_list2)

            df = pd.DataFrame(df_list, columns=headers)

            spa = ' ';

            df['r'] = df['Fraction'].astype(float)/df['Temperature'].astype(float)
            print(df)
            gt=df.sort_values(by=['r'])
            print(gt)
            h = gt['Fraction'].astype(float)
            h2 = gt['Temperature'].astype(float)
            h3 = gt['p'].astype(float)
            h4 = gt['Cp'].astype(float)
            h5 = gt['K'].astype(float)
            h6 = gt['mu'].astype(float)
            h7=  gt['r'].astype(float)

            self.MplWdiget9.canvas.A.scatter(h7, h3, color='black')
            self.MplWdiget9.canvas.A2.scatter(h7, h5, color='black')
            self.MplWdiget9.canvas.A3.scatter(h7, h4, color='black')
            self.MplWdiget9.canvas.A4.scatter(h7, h6, color='black')
            self.MplWdiget9.canvas.draw()

            c1 = np.polyfit(h7, h3, deg=3)
            c2 = np.polyfit(h7, h3, deg=4)
            c3 = np.polyfit(h7, h3, deg=5)

            po = np.poly1d(c1,variable='(phi/T)')
            po2 = np.poly1d(c2,variable='(phi/T)')
            po3 = np.poly1d(c3,variable='(phi/T)')

            s1 = r2_score(h3, po(h7))
            s2 = r2_score(h3, po2(h7))
            s3 = r2_score(h3, po3(h7))

            d1 = np.polyfit(h7, h5, deg=3)
            d2 = np.polyfit(h7, h5, deg=4)
            d3 = np.polyfit(h7, h5, deg=5)

            qo = np.poly1d(d1,variable='(phi/T)')
            qo2 = np.poly1d(d2,variable='(phi/T)')
            qo3 = np.poly1d(d3,variable='(phi/T)')

            t1 = r2_score(h5, qo(h7))
            t2 = r2_score(h5, qo2(h7))
            t3 = r2_score(h5, qo3(h7))

            e1 = np.polyfit(h7, h4, deg=3)
            e2 = np.polyfit(h7, h4, deg=4)
            e3 = np.polyfit(h7, h4, deg=5)

            ro = np.poly1d(e1,variable='(phi/T)')
            ro2 = np.poly1d(e2,variable='(phi/T)')
            ro3 = np.poly1d(e3,variable='(phi/T)')

            u1 = r2_score(h4, ro(h7))
            u2 = r2_score(h4, ro2(h7))
            u3 = r2_score(h4, ro3(h7))

            f1 = np.polyfit(h7, h6, deg=3)
            f2 = np.polyfit(h7, h6, deg=4)
            f3 = np.polyfit(h7, h6, deg=5)

            so = np.poly1d(f1,variable='(phi/T)')
            so2 = np.poly1d(f2,variable='(phi/T)')
            so3 = np.poly1d(f3,variable='(phi/T)')

            w1 = r2_score(h6, so(h7))
            w2 = r2_score(h6, so2(h7))
            w3 = r2_score(h6, so3(h7))



            self.MplWdiget9.canvas.A.plot(h7, po(h7), color='blue')
            self.MplWdiget9.canvas.A.plot(h7, po2(h7), color='red')
            self.MplWdiget9.canvas.A.plot(h7, po3(h7), color='green')

            self.MplWdiget9.canvas.A2.plot(h7, qo(h7), color='blue')
            self.MplWdiget9.canvas.A2.plot(h7, qo2(h7), color='red')
            self.MplWdiget9.canvas.A2.plot(h7, qo3(h7), color='green')

            self.MplWdiget9.canvas.A3.plot(h7, ro(h7), color='blue')
            self.MplWdiget9.canvas.A3.plot(h7, ro2(h7), color='red')
            self.MplWdiget9.canvas.A3.plot(h7, ro3(h7), color='green')

            self.MplWdiget9.canvas.A4.plot(h7, so(h7), color='blue')
            self.MplWdiget9.canvas.A4.plot(h7, so2(h7), color='red')
            self.MplWdiget9.canvas.A4.plot(h7, so3(h7), color='green')
            self.MplWdiget9.canvas.draw()
            self.log_3.setText("")
            self.log_3.append('|###############################Density Correlation################################|')
            self.log_3.append('|----Degreee----|---------r^2--------|')
            self.log_3.append('|1' + 20 * spa + '    |' + str(s1) + '' + 20 * spa + '|')
            self.log_3.append('---------------------------------------')
            self.log_3.append('' + str(po) + '')
            self.log_3.append('---------------------------------------')
            self.log_3.append('|2' + 20 * spa + '     |' + str(s2) + '' + 20 * spa + '|')
            self.log_3.append('---------------------------------------')
            self.log_3.append('' + str(po2) + '')
            self.log_3.append('---------------------------------------')
            self.log_3.append('|3' + 20 * spa + '     |' + str(s3) + '' + 20 * spa + '|')
            self.log_3.append('---------------------------------------')
            self.log_3.append('' + str(po3) + '')
            self.log_3.append('---------------------------------------')
            self.log_3.append('|-----------------------------------------------------------------------------------|')
            self.log_3.append('|########################Thermal Conductivity Correlation###########################|')
            self.log_3.append('|----Degreee----|---------r^2--------|')
            self.log_3.append('|1' + 20 * spa + '|' + str(t1) + '' + 20 * spa + '|')
            self.log_3.append('---------------------------------------')
            self.log_3.append('' + str(ro) + '')
            self.log_3.append('---------------------------------------')
            self.log_3.append('|2' + 20 * spa + '|' + str(t2) + '' + 20 * spa + '|')
            self.log_3.append('---------------------------------------')
            self.log_3.append('' + str(ro2) + '')
            self.log_3.append('---------------------------------------')
            self.log_3.append('|3' + 20 * spa + '|' + str(t3) + '' + 20 * spa + '|')
            self.log_3.append('---------------------------------------')
            self.log_3.append('' + str(ro3) + '')
            self.log_3.append('---------------------------------------')
            self.log_3.append('|-----------------------------------------------------------------------------------|')
            self.log_3.append('|########################Specific Heat Correlation##################################|')
            self.log_3.append('|----Degreee----|---------r^2--------|')
            self.log_3.append('|1' + 20 * spa + '|' + str(u1) + '' + 20 * spa + '|')
            self.log_3.append('---------------------------------------')
            self.log_3.append('' + str(qo) + '')
            self.log_3.append('---------------------------------------')
            self.log_3.append('|2' + 20 * spa + '|' + str(u2) + '' + 20 * spa + '|')
            self.log_3.append('---------------------------------------')
            self.log_3.append('' + str(qo2) + '')
            self.log_3.append('---------------------------------------')
            self.log_3.append('|3' + 20 * spa + '|' + str(u3) + '' + 20 * spa + '|')
            self.log_3.append('---------------------------------------')
            self.log_3.append('' + str(qo3) + '')
            self.log_3.append('---------------------------------------')
            self.log_3.append('|-----------------------------------------------------------------------------------|')
            self.log_3.append('|########################Dynamic Viscocity Correlation##############################|')
            self.log_3.append('|----Degreee----|---------r^2--------|')
            self.log_3.append('|1' + 20 * spa + '|' + str(w1) + '' + 20 * spa + '|')
            self.log_3.append('---------------------------------------')
            self.log_3.append('' + str(so) + '')
            self.log_3.append('---------------------------------------')
            self.log_3.append('|2' + 20 * spa + '|' + str(w2) + '' + 20 * spa + '|')
            self.log_3.append('---------------------------------------')
            self.log_3.append('' + str(so2) + '')
            self.log_3.append('---------------------------------------')
            self.log_3.append('|3' + 20 * spa + '|' + str(w3) + '' + 20 * spa + '|')
            self.log_3.append('---------------------------------------')
            self.log_3.append('' + str(so3) + '')
            self.log_3.append('---------------------------------------')
            self.log_3.append('|-----------------------------------------------------------------------------------|')
        except ValueError:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Warning)
            msg.setWindowTitle("Warning")
            msg.setText("Input Format must be a Numeric Values")

            def msgButtonClick(i):
                msg.close()

            msg.setStandardButtons(QMessageBox.Cancel)
            msg.show()
            msg.buttonClicked.connect(msgButtonClick)

    def adr(self):
        model3 = self.table5.model()
        cxv = self.spin_2.value()
        self.table5.setRowCount(cxv)
    def hyp(self):
        self.tabWidget.setTabVisible(3, True)
        self.tabWidget.setCurrentIndex(3)
        self.spin_2.setValue(0)
        self.log_3.setText("")
        self.table5.setRowCount(0)

        self.MplWdiget9.canvas.A.clear()
        self.MplWdiget9.canvas.A2.clear()
        self.MplWdiget9.canvas.A3.clear()
        self.MplWdiget9.canvas.A4.clear()

        self.MplWdiget9.canvas.A.set_xlabel('phi/T')
        self.MplWdiget9.canvas.A.set_ylabel('Density')
        self.MplWdiget9.canvas.A2.set_xlabel('phi/T')
        self.MplWdiget9.canvas.A2.set_ylabel('Thermal Conductivity')
        self.MplWdiget9.canvas.A3.set_xlabel('phi/T')
        self.MplWdiget9.canvas.A3.set_ylabel('Specific Heat')
        self.MplWdiget9.canvas.A4.set_xlabel('phi/T')
        self.MplWdiget9.canvas.A4.set_ylabel('Dynamic Viscocity')

        self.MplWdiget9.canvas.A.grid(True)
        self.MplWdiget9.canvas.A2.grid(True)
        self.MplWdiget9.canvas.A3.grid(True)
        self.MplWdiget9.canvas.A4.grid(True)
        self.MplWdiget9.canvas.draw()


    def hx(self):
        self.sai.setVisible(False)
    def sx(self):
        self.sai.setVisible(True)

    def hidbox(self):
        self.wok.setVisible(False)

    def showbox(self):
        if self.z3.isChecked()==True:
         self.bill.setVisible(True)
         self.button6.setVisible(False)
         self.wok.setVisible(True)
        else:
         self.wok.setVisible(False)
         self.button6.setVisible(True)

    def abo(self):
        self.tabWidget.setTabVisible(0, True)
        self.tabWidget.setCurrentIndex(0)

    def cle(self):
        self.MplWdiget8.canvas.Axe.set_xlabel('Fraction (%)')
        self.MplWdiget8.canvas.Axe.set_ylabel('Density')
        self.MplWdiget8.canvas.Axe2.set_xlabel('Fraction(%)')
        self.MplWdiget8.canvas.Axe2.set_ylabel('Thermal Conductivity')
        self.MplWdiget8.canvas.Axe3.set_xlabel('Fraction(%)')
        self.MplWdiget8.canvas.Axe3.set_ylabel('Specific Heat')
        self.MplWdiget8.canvas.Axe4.set_xlabel('Fraction(%)')
        self.MplWdiget8.canvas.Axe4.set_ylabel('Dynamic Viscocity')
        self.MplWdiget8.canvas.Axe.grid(True)
        self.MplWdiget8.canvas.Axe2.grid(True)
        self.MplWdiget8.canvas.Axe3.grid(True)
        self.MplWdiget8.canvas.Axe4.grid(True)

        self.MplWdiget8.canvas.Axe.clear()
        self.MplWdiget8.canvas.Axe2.clear()
        self.MplWdiget8.canvas.Axe3.clear()
        self.MplWdiget8.canvas.Axe4.clear()
        self.MplWdiget8.canvas.draw()

    def ro(self):
       try:
        self.MplWdiget8.canvas.Axe.clear()
        self.MplWdiget8.canvas.Axe2.clear()
        self.MplWdiget8.canvas.Axe3.clear()
        self.MplWdiget8.canvas.Axe4.clear()
        self.MplWdiget8.canvas.draw()
        col_count = self.table2.columnCount()
        row_count = self.table2.rowCount()
        headers = [str(self.table2.horizontalHeaderItem(i).text()) for i in range(col_count)]
        df_list = []
        for row in range(row_count):
            df_list2 = []
            for col in range(col_count):
                table_item = self.table2.item(row, col)
                df_list2.append('' if table_item is None else str(table_item.text()))
            df_list.append(df_list2)

        df = pd.DataFrame(df_list, columns=headers)
        spa = ' ';
        x=df['Fraction'].astype(float)
        x2=df['p'].astype(float)
        x3 = df['K'].astype(float)
        x4 = df['Cp'].astype(float)
        x5 = df['mu'].astype(float)

        mo = np.polyfit(x, x2, deg=1)
        mo2= np.polyfit(x,x2,deg=2)
        mo3 = np.polyfit(x, x2, deg=3)

        za=np.polyfit(x, x3, deg=1)
        za2 = np.polyfit(x, x3, deg=2)
        za3 = np.polyfit(x, x3, deg=3)

        sa = np.polyfit(x, x4, deg=1)
        sa2 = np.polyfit(x, x4, deg=2)
        sa3 = np.polyfit(x, x4, deg=3)

        ka = np.polyfit(x, x5, deg=1)
        ka2 = np.polyfit(x, x5, deg=2)
        ka3 = np.polyfit(x, x5, deg=3)

        p = np.poly1d(ka)
        p2 = np.poly1d(ka2)
        p3 = np.poly1d(ka3)

        pr = np.poly1d(sa)
        pr2 = np.poly1d(sa2)
        pr3 = np.poly1d(sa3)

        pri = np.poly1d(za)
        pri2 = np.poly1d(za2)
        pri3 = np.poly1d(za3)

        prid  = np.poly1d(mo)
        prid2 = np.poly1d(mo2)
        prid3 = np.poly1d(mo3)

        sd=   r2_score(x2, prid(x)).__round__(3)
        sd1 = r2_score(x2, prid2(x)).__round__(3)
        sd2 = r2_score(x2, prid3(x)).__round__(3)

        sdq = r2_score(x3, pri(x)).__round__(3)
        sdq1 = r2_score(x3, pri2(x)).__round__(3)
        sdq2 = r2_score(x3, pri3(x)).__round__(3)

        sq = r2_score(x4, pr(x)).__round__(3)
        sq1 = r2_score(x4, pr2(x)).__round__(3)
        sq2 = r2_score(x4, pr3(x)).__round__(3)

        q = r2_score(x5, p(x)).__round__(3)
        q1 = r2_score(x5, p2(x)).__round__(3)
        q2 = r2_score(x5, p3(x)).__round__(3)

        self.log.append('|###############################Density Correlation################################|')
        self.log.append('|----Degreee----|---------r^2--------|')
        self.log.append('|1'+20*spa+'    |'+str(sd)+''+20*spa+'|')
        self.log.append('---------------------------------------')
        self.log.append(''+str(prid)+'')
        self.log.append('---------------------------------------')
        self.log.append('|2'+20*spa+'     |'+str(sd1)+''+20*spa+'|')
        self.log.append('---------------------------------------')
        self.log.append('' + str(prid2) + '')
        self.log.append('---------------------------------------')
        self.log.append('|3'+20*spa+'     |'+str(sd2)+''+20*spa+'|')
        self.log.append('---------------------------------------')
        self.log.append('' + str(prid3) + '')
        self.log.append('---------------------------------------')
        self.log.append('|-----------------------------------------------------------------------------------|')
        self.log.append('|########################Thermal Conductivity Correlation###########################|')
        self.log.append('|----Degreee----|---------r^2--------|')
        self.log.append('|1' + 20 * spa +'|'+ str(sdq) + ''+ 20*spa +'|')
        self.log.append('---------------------------------------')
        self.log.append('' + str(pri) + '')
        self.log.append('---------------------------------------')
        self.log.append('|2' + 20 * spa +'|'+ str(sdq1)+'' +20*spa+'|')
        self.log.append('---------------------------------------')
        self.log.append('' + str(pri2) + '')
        self.log.append('---------------------------------------')
        self.log.append('|3' + 20 * spa +'|'+str(sdq2)+ '' +20*spa+'|')
        self.log.append('---------------------------------------')
        self.log.append('' + str(pri3)+ '')
        self.log.append('---------------------------------------')
        self.log.append('|-----------------------------------------------------------------------------------|')
        self.log.append('|########################Specific Heat Correlation##################################|')
        self.log.append('|----Degreee----|---------r^2--------|')
        self.log.append('|1' + 20 * spa + '|'+ str(sq) + ''+20*spa+'|')
        self.log.append('---------------------------------------')
        self.log.append(''+str(pr)+'')
        self.log.append('---------------------------------------')
        self.log.append('|2'+20*spa+'|'+str(sq1)+''+ 20*spa+'|')
        self.log.append('---------------------------------------')
        self.log.append('' + str(pr2) +'')
        self.log.append('---------------------------------------')
        self.log.append('|3' + 20 * spa +'|'+str(sq2) +''+20*spa+'|')
        self.log.append('---------------------------------------')
        self.log.append(''+str(pr3) +'')
        self.log.append('---------------------------------------')
        self.log.append('|-----------------------------------------------------------------------------------|')
        self.log.append('|########################Dynamic Viscocity Correlation##############################|')
        self.log.append('|----Degreee----|---------r^2--------|')
        self.log.append('|1' + 20 * spa + '|' + str(q) + '' + 20 * spa +'|')
        self.log.append('---------------------------------------')
        self.log.append('' + str(p) + '')
        self.log.append('---------------------------------------')
        self.log.append('|2' + 20 * spa + '|' + str(q1) + '' + 20 * spa +'|')
        self.log.append('---------------------------------------')
        self.log.append('' + str(p2) + '')
        self.log.append('---------------------------------------')
        self.log.append('|3' + 20 * spa + '|' + str(q2) + '' + 20 * spa +'|')
        self.log.append('---------------------------------------')
        self.log.append('' + str(p3) + '')
        self.log.append('---------------------------------------')
        self.log.append('|-----------------------------------------------------------------------------------|')

        self.MplWdiget8.canvas.Axe.set_xlabel('Fraction (%)')
        self.MplWdiget8.canvas.Axe.set_ylabel('Density')
        self.MplWdiget8.canvas.Axe2.set_xlabel('Fraction(%)')
        self.MplWdiget8.canvas.Axe2.set_ylabel('Thermal Conductivity')
        self.MplWdiget8.canvas.Axe3.set_xlabel('Fraction(%)')
        self.MplWdiget8.canvas.Axe3.set_ylabel('Specific Heat')
        self.MplWdiget8.canvas.Axe4.set_xlabel('Fraction(%)')
        self.MplWdiget8.canvas.Axe4.set_ylabel('Dynamic Viscocity')




        self.MplWdiget8.canvas.Axe.scatter(x,x2,color='black',label='d')
        self.MplWdiget8.canvas.Axe.plot(x, prid(x),color='red',label='1')
        self.MplWdiget8.canvas.Axe.plot(x, prid2(x), color='green',label='2')
        self.MplWdiget8.canvas.Axe.plot(x, prid3(x), color='blue',label='3')

        self.MplWdiget8.canvas.Axe2.scatter(x, x3, color='black', label='d')
        self.MplWdiget8.canvas.Axe2.plot(x, pri(x), color='red', label='1')
        self.MplWdiget8.canvas.Axe2.plot(x, pri2(x), color='green', label='2')
        self.MplWdiget8.canvas.Axe2.plot(x, pri3(x), color='blue', label='3')


        self.MplWdiget8.canvas.Axe3.scatter(x, x4, color='black', label='d')
        self.MplWdiget8.canvas.Axe3.plot(x, pr(x), color='red', label='1')
        self.MplWdiget8.canvas.Axe3.plot(x, pr2(x), color='green', label='2')
        self.MplWdiget8.canvas.Axe3.plot(x, pr3(x), color='blue', label='3')


        self.MplWdiget8.canvas.Axe4.scatter(x, x5, color='black', label='d')
        self.MplWdiget8.canvas.Axe4.plot(x, p(x), color='red', label='1')
        self.MplWdiget8.canvas.Axe4.plot(x, p2(x), color='green', label='2')
        self.MplWdiget8.canvas.Axe4.plot(x, p3(x), color='blue', label='3')
        self.MplWdiget8.canvas.Axe.grid(True)
        self.MplWdiget8.canvas.Axe2.grid(True)
        self.MplWdiget8.canvas.Axe3.grid(True)
        self.MplWdiget8.canvas.Axe4.grid(True)
        self.MplWdiget8.canvas.draw()

       except ValueError:
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Warning)
        msg.setWindowTitle("Warning")
        msg.setText("Input Format must be a Numeric Values")

        def msgButtonClick(i):
            msg.close()

        msg.setStandardButtons(QMessageBox.Cancel)
        msg.show()
        msg.buttonClicked.connect(msgButtonClick)

    def rowadd(self):
        model = self.table2.model()

        cx = self.spin.value()
        self.table2.setRowCount(cx)

    def exp(self):

        self.MplWdiget8.canvas.Axe.clear()
        self.MplWdiget8.canvas.Axe2.clear()
        self.MplWdiget8.canvas.Axe3.clear()
        self.MplWdiget8.canvas.Axe4.clear()

        self.spin.setValue(0)
        self.table2.setRowCount(0)
        self.tabWidget.setTabVisible(2, True)
        self.tabWidget.setCurrentIndex(2)

        self.MplWdiget8.canvas.Axe.set_xlabel('Fraction (%)')
        self.MplWdiget8.canvas.Axe.set_ylabel('Density')
        self.MplWdiget8.canvas.Axe2.set_xlabel('Fraction(%)')
        self.MplWdiget8.canvas.Axe2.set_ylabel('Thermal Conductivity')
        self.MplWdiget8.canvas.Axe3.set_xlabel('Fraction(%)')
        self.MplWdiget8.canvas.Axe3.set_ylabel('Specific Heat')
        self.MplWdiget8.canvas.Axe4.set_xlabel('Fraction(%)')
        self.MplWdiget8.canvas.Axe4.set_ylabel('Dynamic Viscocity')
        self.MplWdiget8.canvas.Axe.grid(True)
        self.MplWdiget8.canvas.Axe2.grid(True)
        self.MplWdiget8.canvas.Axe3.grid(True)
        self.MplWdiget8.canvas.Axe4.grid(True)


        self.MplWdiget8.canvas.draw()
        self.log.setText("")

    def ted(self):
      try:
          phi = float(self.sika2.text())
          tem = float(self.v8.text())
          co = np.linspace(0.0, phi / 100, num=11)
          temp = np.linspace(0.0, tem, num=11)
          t = (self.nams.text())
          c = str(t)
          self.MplWidget2.canvas.figure.suptitle("Thermo Physical Properties of " + c)
          model = self.table.model()

          d_w = -4.745 * 1e-6 * pow(temp, 2) - 1.640 * 1e-5 * temp + 1
          k_w = 0.56112 + 0.00193 * temp - 2.601527 * 1e-6 * pow(temp, 2) - 6.08803 * 1e-8 * pow(temp, 3)
          mu_w = 0.00169 - 4.252 * 1e-5 * temp + 4.9255 * 1e-7 * pow(temp, 2) - 2.09935 * 1e-9 * pow(temp, 3)
          c_w = 4217.629 - 3.2088 * temp + 0.09503 * pow(temp, 2) - 0.00132 * pow(temp, 3) + 9.415 * 1e-6 * pow(temp,4) - 2.5479 * 1e-8 * pow(temp, 5)

          val = model.index(0, 0)
          val1 = model.index(0, 1)
          val2 = model.index(0, 2)
          val3 = model.index(0, 3)

          va = model.index(1, 0)
          va1 = model.index(1, 1)
          va2 = model.index(1, 2)
          va3 = model.index(1, 3)

          v = model.index(2, 0)
          v1 = model.index(2, 1)
          v2 = model.index(2, 2)
          v3 = model.index(2, 3)


          p1 = float(model.data(val))
          p2 = float(model.data(val1))
          p3 = float(model.data(val2))
          p4 = float(model.data(val3))

          c1 = float(model.data(va))
          c2 = float(model.data(va1))
          c3 = float(model.data(va2))
          c4 = float(model.data(va3))

          b1 = float(model.data(v))
          b2 = float(model.data(v1))
          b3 = float(model.data(v2))
          b4 = float(model.data(v3))


          d_n = p1 + (p2) * temp + (p3) * pow(temp, 2) + (p4) * pow(temp, 3)
          c_p = c1 + (c2) * temp + (c3) * pow(temp, 2) + (c4) * pow(temp, 3)
          k_p = b1 + (b2) * temp + (b3) * pow(temp, 2) + (b4) * pow(temp, 3)
          x1 = float(self.v5.text())
          x2 = float(self.v6.text())
          x3 = float(self.v7.text())

          d_nf = (1 - co) * d_w + co * d_n
          cp_nf = (1 - co) * c_w + co * c_p
          k_nf = ((k_p + 2 * k_w + 2 * co * (k_p - k_w)) / (k_p + 2 * k_w - co * (k_p - k_w))) * k_w
          mu_nf = mu_w * (x1 + x2 * phi + x3 * pow(phi, 2))

          self.MplWidget2.canvas.Axes3D.plot(temp, co, d_nf, 'darkorange')
          self.MplWidget2.canvas.Axes3D2.plot(temp, co, cp_nf, 'darkgreen')
          self.MplWidget2.canvas.Axes3D3.plot(temp, co, k_nf, 'darkblue')
          self.MplWidget2.canvas.Axes3D4.plot(temp, co, mu_nf, 'darkviolet')
          self.MplWidget2.canvas.draw()
      except ValueError:
           msg = QMessageBox()
           msg.setIcon(QMessageBox.Warning)
           msg.setWindowTitle("Warning")
           msg.setText("Chek the input format")
           def msgButtonClick(i):
             msg.close()

           msg.setStandardButtons(QMessageBox.Cancel)
           msg.show()
           msg.buttonClicked.connect(msgButtonClick)




    def texs(self):
        self.tabWidget.setTabVisible(4,True)
        self.tabWidget.setCurrentIndex(4)
        self.nams.setText("")
        self.sika2.setText("")
        self.v8.setText("")
        self.wok.setVisible(False)
        self.z3.setChecked(False)

        self.MplWidget2.canvas.Axes3D.clear()
        self.MplWidget2.canvas.Axes3D2.clear()
        self.MplWidget2.canvas.Axes3D3.clear()
        self.MplWidget2.canvas.Axes3D4.clear()

        self.MplWidget2.canvas.Axes3D.grid(True)
        self.MplWidget2.canvas.Axes3D2.grid(True)
        self.MplWidget2.canvas.Axes3D3.grid(True)
        self.MplWidget2.canvas.Axes3D4.grid(True)

        self.MplWidget2.canvas.Axes3D.set_xlabel("Temperature (C)")
        self.MplWidget2.canvas.Axes3D.set_ylabel("Concentration (%)")
        self.MplWidget2.canvas.Axes3D.set_zlabel("Density (Kg/m^3)")

        self.MplWidget2.canvas.Axes3D2.set_xlabel("Temperature (C)")
        self.MplWidget2.canvas.Axes3D2.set_ylabel("Concentration (%)")
        self.MplWidget2.canvas.Axes3D2.set_zlabel("Specific heat (Kg/m.K)")

        self.MplWidget2.canvas.Axes3D3.set_xlabel("Temperature (C)")
        self.MplWidget2.canvas.Axes3D3.set_ylabel("Concentration (%)")
        self.MplWidget2.canvas.Axes3D3.set_zlabel("Thermal Conductivity (W/m.K)")

        self.MplWidget2.canvas.Axes3D4.set_xlabel("Temperature (C)")
        self.MplWidget2.canvas.Axes3D4.set_ylabel("Concentration (%)")
        self.MplWidget2.canvas.Axes3D4.set_zlabel("Viscoscity (Kg/m.s)")
        self.MplWidget2.canvas.draw()


    def act(self):
       try:
        te = (self.a1.text())
        con = str(te)
        self.MplWidget.canvas.figure.suptitle("Physical Properties of " + con)
        phi = float(self.v1.text())
        den = float(self.v2.text())
        cp = float(self.v3.text())
        k_p = float(self.v4.text())
        x1 = float(self.v5.text())
        x2 = float(self.v6.text())
        x3 = float(self.v7.text())
        p_w = 998
        k_w = 0.643
        c_w = 4189.0
        mu_w = 0.00089
        if self.t1.isChecked() == True:
            p_w = 998
            k_w = 0.643
            c_w = 4189.0
            mu_w = 1.5e-3
        if self.t2.isChecked()==True:
            p_w = 1110
            k_w = 0.258
            c_w = 3583
            mu_w = 0.0162
        if self.t3.isChecked()==True:
            p_w = float(self.ba.text())
            k_w = float(self.ba3.text())
            c_w = float(self.ba2.text())
            mu_w = float(self.ba4.text())

        X_arr = np.linspace(0.0, phi / 100, num=11)
        y1 = (1 - X_arr) * p_w + X_arr * den
        self.MplWidget.canvas.axes.plot(X_arr * 100, y1, color='blue')
        y2 = (1 - X_arr) * c_w + X_arr * cp
        self.MplWidget.canvas.axes2.plot(X_arr * 100, y2, color='red')
        y3 =  mu_w* (1+x1*X_arr + x2 * X_arr*X_arr + x3 * X_arr * X_arr*X_arr)
        self.MplWidget.canvas.axes4.plot(X_arr * 100, y3, color='green')
        y4 = ((k_p + 2 * k_w + 2 * X_arr * (k_p - k_w)) / (k_p + 2 * k_w - X_arr * (k_p - k_w))) * k_w
        self.MplWidget.canvas.axes3.plot(X_arr * 100, y4, color='magenta')
        self.MplWidget.canvas.draw()
        self.jup.setEnabled(False)
       except ValueError:
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Warning)
        msg.setWindowTitle("Warning")
        msg.setText("Input Format is incorrect")
        def msgButtonClick(i):
         msg.close()
         msg.setStandardButtons(QMessageBox.Cancel)
         msg.show()
         msg.buttonClicked.connect(msgButtonClick)

    def viz(self):
        self.MplWidget.canvas.axes.clear()
        self.MplWidget.canvas.axes2.clear()
        self.MplWidget.canvas.axes3.clear()
        self.MplWidget.canvas.axes4.clear()
        self.MplWidget.canvas.axes.set_xlabel("Concentration (%)")
        self.MplWidget.canvas.axes.set_ylabel("Density (Kg/m^3)")

        self.MplWidget.canvas.axes2.set_xlabel("Concentration (%)")
        self.MplWidget.canvas.axes2.set_ylabel("Specific heat (Kg/m.K)")

        self.MplWidget.canvas.axes3.set_xlabel("Concentration (%)")
        self.MplWidget.canvas.axes3.set_ylabel("Thermal Conductivity (W/m.K)")

        self.MplWidget.canvas.axes4.set_xlabel("Concentration (%)")
        self.MplWidget.canvas.axes4.set_ylabel("Viscoscity (Kg/m.s)")

        self.MplWidget.canvas.axes.grid(True)
        self.MplWidget.canvas.axes2.grid(True)
        self.MplWidget.canvas.axes3.grid(True)
        self.MplWidget.canvas.axes4.grid(True)
        self.MplWidget.canvas.draw()

        self.tabWidget.setCurrentIndex(1)
        self.tabWidget.setTabVisible(1, True)
        self.jup.setEnabled(True)

    def print(self):
       try:
        self.tex.setText("")
        self.jo2.setChecked(False)
        self.jbox.setVisible(False)
        self.tabWidget.setCurrentIndex(5)
        self.tabWidget.setTabVisible(5, True)
        self.button3.setEnabled(True)
        self.button3.setVisible(True)

        te = (self.a1.text())
        con = str(te)
        k_p = float(self.v4.text())
        phi = float(self.v1.text())
        den = float(self.v2.text())
        cp = float(self.v3.text())
        k_p = float(self.v4.text())
        x1 = float(self.v5.text())
        x2 = float(self.v6.text())
        x3 = float(self.v7.text())
        k_w = 0.643
        x = ' '
        self.tex.append('/****************ANSYS FLUENT UDF FOR '+(con)+ +1*x+'NANOFLUID****************/')
        self.tex.append('/****Compile before starting the Fluent****/')
        self.tex.append('#include "udf.h"')
        self.tex.append('#define k_p %f '%(k_p))
        self.tex.append('#define phi %f ' % (phi/100))
        self.tex.append('#define cpp %f ' % (cp))
        self.tex.append('#define rho_p %f ' % (den))
        self.tex.append('DEFINE_PROPERTY(nanofluid_conductivity,c,t)')
        self.tex.append('{')
        self.tex.append('real k_w;')
        self.tex.append('real k_nf;')
        self.tex.append('real ctemp = C_T(c,t);')
        self.tex.append('k_w=-0.000008*(ctemp*ctemp)+((0.0062*ctemp)-0.5388);')
        self.tex.append('k_nf=(((k_p+(2*k_w)-(2*(k_w-k_p))*phi)/(k_p+(2*k_w)+(k_w-k_p)*phi)))*k_w;')
        self.tex.append('return k_nf;')
        self.tex.append('}')
        self.tex.append('DEFINE_PROPERTY(nanofluid_viscosity,c,t)')
        self.tex.append('{')
        self.tex.append('real mu_w;')
        self.tex.append('real mu_nf;')
        self.tex.append('real ctemp2 = C_T(c,t);')
        self.tex.append('mu_w = 2.414*0.00001*pow(10,(247.8/(ctemp2-140)));')
        self.tex.append('mu_nf = mu_w*((%f)+(%f)*phi+(%f)*pow(phi,2));'%(x1,x2,x3))
        self.tex.append('return mu_nf;')
        self.tex.append('}')
        self.tex.append('DEFINE_SPECIFIC_HEAT(nanofluid_specificheat,c,t)')
        self.tex.append('{')
        self.tex.append('real cp_w;')
        self.tex.append('real cp_nf;')
        self.tex.append('real ctemp3 = C_T(c,t);')
        self.tex.append('cp_w=4.21-0.00561*ctemp3+0.001299*pow(ctemp3,1.5)-0.0001153*pow(ctemp3,2)+0.0006965*pow(ctemp3,2.5);')
        self.tex.append('cp_nf = (1-phi)*cp_w+phi*cpp;')
        self.tex.append('return cp_nf;')
        self.tex.append('}')
        self.tex.append('DEFINE_PROPERTY(nanofluid_density,c,t)')
        self.tex.append('{')
        self.tex.append('real rho_w;')
        self.tex.append('real rho_nf;')
        self.tex.append('real ctemp4 = C_T(c,t);')
        self.tex.append('rho_w = (-3.570*(pow(10,-3))*(pow(ctemp4,2))+(1.88*ctemp4+753));')
        self.tex.append('rho_nf = (1-phi)*rho_w+phi*rho_p;')
        self.tex.append('}')

       except ValueError:
         self.button3.setEnabled(False)
         msg = QMessageBox()
         msg.setIcon(QMessageBox.Warning)
         msg.setWindowTitle("Warning")
         msg.setText("UDF cannot be save ")
         def msgButtonClick(i):
          msg.close()
         msg.setStandardButtons(QMessageBox.Cancel)
         msg.show()
         msg.buttonClicked.connect(msgButtonClick)

    def save(self):
        filename=QFileDialog.getSaveFileName(self,'Save File',os.getenv('HOME'))
        with open (filename[0],'w') as f:
            text=self.tex.toPlainText()
            f.write(text)
            self.tex.append("------------------------------")
            self.tex.append("File Saved")
            self.tex.append("------------------------------")

if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    app.setStyle('fusion')
    window = MyWindow()
    window.show()
    sys.exit(app.exec_())