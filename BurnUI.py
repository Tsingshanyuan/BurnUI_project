# -*- coding: utf-8 -*-
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from PyQt5.QtWidgets import QVBoxLayout, QMainWindow, QWidget, QLabel, QComboBox, \
    QTabWidget, QMessageBox, QAction, QTextBrowser, QFileDialog, QPushButton,  QHBoxLayout, \
    QLineEdit, QApplication, QCheckBox, QDialog, QGridLayout, QProgressBar
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QIcon
from scipy.integrate import odeint
from numpy import array, append, interp, arange, linspace
from matplotlib import pyplot as plt
from sys import argv


class Nuclide():
    """表示某个核素的类"""

    def __init__(self, name, sig_a, sig_g, parent=0, sig_f=0, T=-1, lamda=0.0, N0=0):
        """初始化"""
        self.name = name
        self.sig_a = sig_a*barn
        self.sig_a_origin = self.sig_a
        self.sig_a_input = QLineEdit(str(self.sig_a/barn))
        self.sig_a_input.setToolTip('原始数据(0.0253eV)：' + str(sig_a))
        self.sig_g = sig_g*barn
        self.sig_g_origin = self.sig_g
        self.sig_g_input = QLineEdit(str(self.sig_g/barn))
        self.sig_g_input.setToolTip('原始数据(0.0253eV)：' + str(sig_g))
        self.sig_f = sig_f*barn
        self.sig_f_origin = self.sig_f
        self.sig_f_input = QLineEdit(str(self.sig_f/barn))
        self.sig_f_input.setToolTip('原始数据(0.0253eV)：' + str(sig_f))
        if lamda == 0 and T == -1:
            self.lamda = 0
        elif lamda == 0:
            self.lamda = 0.6931471805599453/T
        elif T == -1:
            self.lamda = lamda
        self.lamda_origin = self.lamda
        self.lamda_input = QLineEdit(str('%.4e' % self.lamda))
        self.lamda_input.setToolTip('原始数据：' + str(lamda))
        self.N0 = N0*N
        self.Nt = array([self.N0])
        self.parent = parent
        self.Show_Nuc()
        Nuclides[self.name] = self

    def solve_Nt(self, f, t, j1, P, V):
        self.Nt = append(self.Nt, odeint(f, self.Nt[-1], t, args=(self.parent, self, P, V, j1))[-1])

    def Sigma_f(self, j1):
        return(self.sig_f*self.Nt[j1-1])

    def Nplot(self, t, type):
        if type == 0:
            plt.loglog(t, self.Nt, label=self.name, linewidth=0.5)
        elif type == 1:
            plt.plot(t, self.Nt, label=self.name, linewidth=0.5)

    def get_input(self):
        self.N0 = float(self.input.text())
        self.Nt = array([self.N0])

    def show_output(self, t, ti):
        result = self.Ninterp(t, ti)
        self.output.setText(str('%.4e' % result))

    def Ninterp(self, t, ti):
        return(interp(t, ti, self.Nt))

    def Show_Nuc(self):
        self.label = QCheckBox('---' + self.name + ":")
        self.label.setToolTip("选择是否显示曲线")
        self.label.setChecked(True)
        self.input = QLineEdit(str('%.4e' % self.N0))
        self.input.setToolTip("请输入该核素的初始核密度/10^22cm^-3")
        self.output = QLineEdit("0.0")


class Product(Nuclide):
    def __init__(self, name, sig_a, sig_g, gama25, gama49, gama41, parent=0, sig_f=0, T=-1, lamda=0, N0=0):
        super().__init__(name, sig_a, sig_g, parent=parent, sig_f=sig_f, T=T, lamda=lamda, N0=N0)

        # 裂变产额
        self.gama25 = gama25/100
        self.gama25_origin = gama25/100
        self.gama25_input = QLineEdit(str(gama25))
        self.gama25_input.setToolTip('原始数据：' + str(gama25))
        self.gama49 = gama49/100
        self.gama49_origin = gama49 / 100
        self.gama49_input = QLineEdit(str(gama49))
        self.gama49_input.setToolTip('原始数据：' + str(gama49))
        self.gama41 = gama41/100
        self.gama41_origin = gama41 / 100
        self.gama41_input = QLineEdit(str(gama41))
        self.gama41_input.setToolTip('原始数据：' + str(gama41))


def getRealPath(s):
 # 获取exe解压目录的绝对路径
 from os.path import realpath
 from sys import path
 p = realpath(path[0])
 p = p.replace(r'\base_library.zip', '')
 p = p + s
 return p


def f(y, t, U1, U2, P, V, j1):
    # 燃耗方程，根据模式而变
    if s == 0:  # 定通量模式
        phi = P
    elif s == 1:  # 定功率模式
        Sigma = (Nuclides["U235 "].Sigma_f(j1) + Nuclides["Pu239"].Sigma_f(j1) + Nuclides["Pu241"].Sigma_f(j1))
        phi = 3.12e10 * P / (V * Sigma)
    if type(U2) == type(I135):
        sigma_phi = (U2.gama25 * Nuclides["U235 "].Sigma_f(j1) + U2.gama49 * Nuclides["Pu239"].Sigma_f(j1) +
                      U2.gama41 * Nuclides["Pu241"].Sigma_f(j1))*phi
    elif type(U2) == type(U235):
        sigma_phi = U1.sig_g*phi*U1.Nt[j1]
    return array(sigma_phi + U1.lamda*U1.Nt[j1] - (U2.sig_a*phi + U2.lamda)*y)


class Main_Widget(QWidget):
    # 窗口内容
    def __init__(self):
        super().__init__()
        # 总时间
        self.t_tot = array([0])

        # 定义QWidgets
        self.figure = plt.figure(dpi=300)
        self.canvas = FigureCanvas(self.figure)
        self.dataBrowser = QTextBrowser()
        self.dataWidget = QTabWidget()
        self.canvasTab = QWidget()
        self.dataTab = QWidget()
        self.canvasbox = QVBoxLayout(self.canvasTab)
        self.databox = QVBoxLayout(self.dataTab)
        self.canvasbox.addWidget(self.canvas)
        self.databox.addWidget(self.dataBrowser)
        self.dataWidget.addTab(self.canvasTab, "图形")
        self.dataWidget.addTab(self.dataTab, "数据")

        self.button_draw = QPushButton("计算")
        self.button_draw.setShortcut(Qt.Key_Return)
        self.button_draw.setToolTip("计算并绘图(快捷键:回车)")
        self.hbox_mode = QHBoxLayout()
        self.lbl_mode = QLabel("选择模式:")
        self.combo_mode = QComboBox()
        self.combo_mode.addItem("定通量燃耗")
        self.combo_mode.addItem("定功率燃耗")
        self.hbox_mode.addWidget(self.lbl_mode)
        self.hbox_mode.addWidget(self.combo_mode)
        self.ipt_phi, self.hbox_phi = self.add_input("中子通量:  ", "1e13", "/cm^2/s")
        self.ipt_Dt, self.hbox_Dt = self.add_input("外步长:    ", "30", "day    ")
        self.ipt_Dt.setToolTip("指整个燃耗时间段的划分，每点击一次‘计算’，程序自动计算一个外步长的时间，并停下来等待用户修改通量、功率等参数")
        self.ipt_dt, self.hbox_dt = self.add_input("内步长:    ", "1", "day    ")
        self.ipt_dt.setToolTip("指燃耗方程数值求解过程中的步长，影响程序计算的精度")
        self.ipt_P, self.hbox_P = self.add_input("反应堆功率:", "300", "MW     ")
        self.ipt_V, self.hbox_V = self.add_input("反应堆体积:", "16", "m^3    ")
        self.tabWidget1 = QTabWidget()
        self.tabWidget2 = QTabWidget()
        self.tab1 = QWidget()
        self.tab2 = QWidget()
        self.tab3 = QWidget()
        self.tabWidget1.addTab(self.tab1, "初值")
        self.tabWidget1.addTab(self.tab2, "结果")
        self.tabWidget2.addTab(self.tab3, "显示")
        self.lbl_plotType = QLabel("坐标类型:")
        self.hbox_plotType = QHBoxLayout()
        self.combo_plotType = QComboBox()
        self.combo_plotType.addItem("对数")
        self.combo_plotType.addItem("线性")
        self.title_Nuc = QLabel("核密度 cm^-3     ")
        self.hbox_plotType.addWidget(self.lbl_plotType)
        self.hbox_plotType.addWidget(self.combo_plotType)
        self.hbox_plotType.addWidget(self.title_Nuc)
        self.ipt_t_search = QLineEdit(str(self.t_tot[-1]))
        self.ipt_t_search.setToolTip("输入查询的时刻/day")
        self.lbl_t_search = QLabel("t = ")
        self.unit_t_search = QLabel("day")
        self.btn_t_search = QPushButton("查询")
        self.btn_t_search.setToolTip("点击查询t时刻的核密度")
        self.btn_t_search.clicked.connect(self.t_search)
        self.pbar = QProgressBar()
        self.pbar.setStyleSheet("QProgressBar{border:1px solid Silver;"
                                "height:20;"
                                "background: white;"
                                "text-align:center;}"

                                "QProgressBar::chunk{background-color:Yellow;"
                                "width:1px;}"  # 宽度
                                )
        self.pbar.setValue(0)

        # 空行
        self.lbl_blank1 = QLabel(' ')

        # 复选框
        self.combo_mode.currentIndexChanged.connect(self.enable_inputs)
        self.combo_mode.currentIndexChanged.connect(self.change_status)
        self.enable_inputs()

        # 参数初值
        s = self.combo_mode.currentIndex()
        self.s = array([])
        self.Pphi = array([])
        self.V = array([])

        # 连接事件
        self.button_draw.clicked.connect(self.Draw)
        self.combo_plotType.currentIndexChanged.connect(self.show_Nuclides)
        for name, Nuc in Nuclides.items():
            Nuc.label.clicked.connect(self.show_Nuclides)

        # 设置布局
        vbox = QVBoxLayout()
        vbox.addLayout(self.hbox_mode)
        vbox.addLayout(self.hbox_phi)
        vbox.addLayout(self.hbox_P)
        vbox.addLayout(self.hbox_V)
        vbox.addLayout(self.hbox_Dt)
        vbox.addLayout(self.hbox_dt)
        vbox.addWidget(self.lbl_blank1)
        vbox.addLayout(self.hbox_plotType)
        vbox1 = QVBoxLayout()
        vbox2 = QVBoxLayout()
        vbox3 = QVBoxLayout()
        hbox1 = QHBoxLayout()
        hbox2 = QHBoxLayout()
        hbox2.addWidget(self.lbl_t_search)
        hbox2.addWidget(self.ipt_t_search)
        hbox2.addWidget(self.unit_t_search)
        hbox2.addWidget(self.btn_t_search)
        for name, Nuc in Nuclides.items():
            vbox1.addWidget(Nuc.label)
            vbox2.addWidget(Nuc.input)
            vbox3.addWidget(Nuc.output)
        self.tab1.setLayout(vbox2)
        self.tab2.setLayout(vbox3)
        self.tab3.setLayout(vbox1)
        hbox1.addWidget(self.tabWidget2)
        hbox1.addWidget(self.tabWidget1)
        vbox.addLayout(hbox1)
        vbox.addLayout(hbox2)
        vbox.addWidget(self.button_draw)
        vbox.addWidget(self.pbar)
        hbox = QHBoxLayout()
        hbox.addLayout(vbox)
        hbox.addWidget(self.dataWidget)
        self.setLayout(hbox)

    def Draw(self):
        self.button_draw.setDisabled(True)
        self.pbar.setStyleSheet("QProgressBar{border:1px solid Silver;"
                                "height:20;"
                                "background: white;"
                                "text-align:center;}"

                                "QProgressBar::chunk{background-color:Yellow;"
                                "width:1px;}"  # 宽度
                                )
        try:
            # 输入参数
            Dt = float(self.ipt_Dt.text()) * day
            dt = float(self.ipt_dt.text()) * day
            # 输入初始核密度
            if len(self.t_tot) == 1:
                for name, Nuc in Nuclides.items():
                    Nuc.get_input()
            # 计算时间步
            j = int(Dt // dt) + 1
            ti = arange(0, Dt + dt, dt)
            self.t_tot = append(self.t_tot, self.t_tot[-1] + ti[1::])
            self.ipt_t_search.setText(str(self.t_tot[-1] / day))

            s = self.combo_mode.currentIndex()
            if len(self.Pphi) == 0:
                self.V = append(self.V, float(self.ipt_V.text()) * 1e4)
                self.s = append(self.s, s)
                if s == 0:
                    self.Pphi = append(self.Pphi, float(self.ipt_phi.text()))
                elif s == 1:
                    self.Pphi = append(self.Pphi, float(self.ipt_P.text())*1e6)
            self.V = append(self.V, float(self.ipt_V.text()) * 1e4 + 0 * ti[1::])
            self.s = append(self.s, s + 0*ti[1::])
            if s == 0:
                self.Pphi = append(self.Pphi, float(self.ipt_phi.text()) + 0*ti[1::])
            elif s == 1:
                self.Pphi = append(self.Pphi, float(self.ipt_P.text())*1e6 + 0*ti[1::])


            # 计算各核素的核密度随时间的变化
            U0.Nt = 0 * self.t_tot
            tstart = len(self.t_tot) - len(ti)
            for j1 in tstart + arange(1, j):
                self.pbar.setValue(int(100*(j1-tstart)/(j-1)))
                t = linspace(self.t_tot[j1 - 1], self.t_tot[j1], 2)
                for name, Nuc in Nuclides.items():
                    Nuc.solve_Nt(f, t, j1, self.Pphi[j1], self.V[j1])
                    app.processEvents()

            # 输出最终核密度
            for name, Nuc in Nuclides.items():
                Nuc.show_output(self.t_tot[-1], self.t_tot)
                app.processEvents()
            self.show_Nuclides()
        except ValueError:
            QMessageBox.warning(self, '警告', "输入错误，请重新输入！", QMessageBox.Ok, QMessageBox.Ok)
        self.button_draw.setEnabled(True)


    def enable_inputs(self):
        # 选择模式“定通量”、“定功率”
        s = self.combo_mode.currentIndex()
        if s == 0:
            self.ipt_phi.setEnabled(True)
            self.ipt_P.setDisabled(True)
            self.ipt_V.setDisabled(True)
        elif s == 1:
            self.ipt_phi.setDisabled(True)
            self.ipt_P.setEnabled(True)
            self.ipt_V.setEnabled(True)

    def change_status(self):
        s = self.combo_mode.currentIndex()
        if s == 0:
            main_window.statusBar().showMessage('-定通量燃耗-  '+'时刻: ' + str(self.t_tot[-1] / day) + 'day' + '  中子通量: ' +
                                                self.ipt_phi.text() + "/cm^2/s")
        elif s == 1:
            main_window.statusBar().showMessage('-定功率燃耗-  '+'时刻: ' + str(self.t_tot[-1] / day) + 'day' + '  反应堆功率: ' +
                                                self.ipt_P.text() + "MW" + '  反应堆体积: ' + self.ipt_V.text()
                                                + "m^3")

    def t_search(self):
        # 查询特定时间的核密度
        self.tabWidget1.setCurrentIndex(1)
        t = float(self.ipt_t_search.text()) * day
        if t >= 0 and t <= self.t_tot[-1]:
            for name, Nuc in Nuclides.items():
                Nuc.show_output(t, self.t_tot)
        else:
            QMessageBox.warning(self, '警告', "查询超出数据范围！", QMessageBox.Ok, QMessageBox.Ok)

    def add_input(self, label, input, unit):
        hbox = QHBoxLayout()
        lbl = QLabel(label)
        ipt = QLineEdit(input)
        unt = QLabel(unit)
        hbox.addWidget(lbl)
        hbox.addWidget(ipt)
        hbox.addWidget(unt)
        return ipt, hbox

    def show_Nuclides(self):
        self.pbar.setStyleSheet("QProgressBar{border:1px solid Silver;"
                                "height:20;"
                                "background: Yellow;"
                                "text-align:center;}"

                                "QProgressBar::chunk{background-color:OrangeRed;"
                                "width:1px;}"  # 宽度
                                )
        if len(self.t_tot) != 1:
            self.figure.clear()

            self.dataText = '核密度(cm-3)随时间的变化\n'
            self.dataText += '时间(d)'
            self.dataText += '\tphi/cm2s\t功率(MW)\t体积(m3)'
            for name, Nuc in Nuclides.items():
                if Nuc.label.isChecked():
                    self.dataText += '\t' + Nuc.name
                app.processEvents()
            self.dataText += '\n'
            for i1 in arange(0, len(self.t_tot)):
                self.dataText += str(self.t_tot[i1]/day)
                if self.s[i1] == 0:  # 定通量模式
                    showphi = '\t' + str('%.2e' % self.Pphi[i1]) + '\t--\t--'
                elif self.s[i1] == 1:  # 定功率模式
                    showPphi = self.Pphi[i1]/1e6
                    showV = self.V[i1]/1e4
                    showphi = '\t--\t' + str('%.2e' % showPphi) + '\t' + str(showV)
                self.dataText += showphi
                for name, Nuc in Nuclides.items():
                    if Nuc.label.isChecked():
                        self.dataText += '\t' + str('%.2e' % Nuc.Nt[i1])
                    app.processEvents()
                self.dataText += '\n'
                self.pbar.setValue(int(100 * i1 / (len(self.t_tot)-1)))
                app.processEvents()

            for name, Nuc in Nuclides.items():
                if Nuc.label.isChecked():
                    Nuc.Nplot(self.t_tot/day, self.combo_plotType.currentIndex())
                    app.processEvents()
            plt.legend(fontsize=4)
            plt.title('Atomic Densities - Time', fontsize=6)
            plt.xlabel("Time/day", fontsize=6)
            plt.ylabel('Atomic Densities/$cm^-3$', fontsize=6)
            plt.xticks(fontsize=4)
            plt.yticks(fontsize=4)
            # 列表
            self.dataBrowser.setText(self.dataText)
            # 画图
            self.canvas.draw()
            self.tabWidget1.setCurrentIndex(1)
            self.change_status()

    def clear(self):
        reply = QMessageBox.question(self, '消息',"确认清空数据？", QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        if reply == QMessageBox.Yes:
            self.pbar.setStyleSheet("QProgressBar{border:1px solid Silver;"
                                    "height:20;"
                                    "background: white;"
                                    "text-align:center;}"

                                    "QProgressBar::chunk{background-color:Yellow;"
                                    "width:1px;}"  # 宽度
                                    )
            self.pbar.setValue(0)
            self.t_tot = array([0])
            self.Pphi = array([])
            self.V = array([])
            self.s = array([])
            self.tabWidget1.setCurrentIndex(0)
            self.ipt_t_search.setText("0")
            for name, Nuc in Nuclides.items():
                Nuc.get_input()
                Nuc.output.setText("0.0")
            self.figure.clear()
            self.canvas.draw()
            self.dataBrowser.setText('')
        else:
            pass


class Main_window(QMainWindow):
    # 主程序窗口
    def __init__(self):
        super().__init__()

        self.initUI()

    def initUI(self):
        self.Draw_Widget = Main_Widget()
        self.setCentralWidget(self.Draw_Widget)
        self.setWindowTitle("BurnUI-点燃耗计算程序")
        path = '\Firecolor.ico'
        self.setWindowIcon(QIcon(getRealPath(path)))
        s = self.Draw_Widget.combo_mode.currentIndex()
        if s == 0:
            self.statusBar().showMessage('-定通量燃耗-  '+'时刻: '+str(self.Draw_Widget.t_tot[-1]/day)+'day'+ '  中子通量: ' +
                                                self.Draw_Widget.ipt_phi.text() + "/cm^2/s")
        elif s == 1:
            self.statusBar().showMessage('-定功率燃耗-  '+'时刻: ' + str(self.Draw_Widget.t_tot[-1] / day) + 'day' + '  反应堆功率: ' +
                                         self.Draw_Widget.ipt_P.text() + "MW" + '  反应堆体积: ' + self.Draw_Widget.ipt_V.text()
                                         + "m^3")
        
        path = '\Clear.ico'
        clearAction = QAction(QIcon(getRealPath(path)), '清空数据, Ctrl+Shift+C', self)
        clearAction.setShortcut('Ctrl+Shift+C')
        clearAction.triggered.connect(self.Draw_Widget.clear)
        path = '\save.ico'
        savedataAction = QAction(QIcon(getRealPath(path)), '保存数据, Ctrl+S', self)
        savedataAction.setShortcut('Ctrl+S')
        savedataAction.triggered.connect(self.saveData)
        path = '\save2.ico'
        savefigureAction = QAction(QIcon(getRealPath(path)), '保存图形, Ctrl+Shift+S', self)
        savefigureAction.setShortcut('Ctrl+Shift+S')
        savefigureAction.triggered.connect(self.saveFigure)
        path = '\edit.ico'
        editAction = QAction(QIcon(getRealPath(path)), '修改核数据', self)
        editAction.triggered.connect(self.edit)
        path = '\help.ico'
        helpAction = QAction(QIcon(getRealPath(path)), '帮助', self)
        helpAction.triggered.connect(self.help)
        path = '\info.ico'
        infoAction = QAction(QIcon(getRealPath(path)), '关于', self)
        infoAction.triggered.connect(self.show_info)
        self.toolbar = self.addToolBar('Clear')
        self.toolbar.addAction(savedataAction)
        self.toolbar.addAction(savefigureAction)
        self.toolbar.addAction(clearAction)
        self.toolbar.addAction(editAction)
        self.toolbar.addAction(helpAction)
        self.toolbar.addAction(infoAction)

        self.edit_window = QDialog(self)
        self.edit_window.setWindowTitle('修改核数据')
        self.edit_window.setGeometry(700, 200, 500, 300)
        grid = QGridLayout()
        self.edit_window.setLayout(grid)

        # 修改和数据的弹窗
        titles = ['核素', '吸收截面/b', '俘获截面/b', '裂变截面/b', '衰变常数/s-1']
        for i1 in range(5):
            title = QLabel(titles[i1])
            grid.addWidget(title, *(0, i1))
        i = 1
        for name, Nuc in Nuclides.items():
            label = QLabel(name)
            grid.addWidget(label, *(i, 0))
            grid.addWidget(Nuc.sig_a_input, *(i, 1))
            grid.addWidget(Nuc.sig_g_input, *(i, 2))
            grid.addWidget(Nuc.sig_f_input, *(i, 3))
            grid.addWidget(Nuc.lamda_input, *(i, 4))
            i += 1
        title2 = ['裂变产物', 'U235 ', 'Pu239', 'Pu241']
        for i1 in range(4):
            title = QLabel(title2[i1])
            grid.addWidget(title, *(i, i1))
        i += 1
        for name, Nuc in Nuclides.items():
            if type(Nuc) == type(I135):
                label = QLabel(name)
                unit = QLabel('(裂变产额%)')
                grid.addWidget(label, *(i, 0))
                grid.addWidget(Nuc.gama25_input, *(i, 1))
                grid.addWidget(Nuc.gama49_input, *(i, 2))
                grid.addWidget(Nuc.gama41_input, *(i, 3))
                grid.addWidget(unit, *(i, 4))
                i += 1
        note = QLabel("注：原始数据来源于ENDF，能量0.0253eV")
        grid.addWidget(note, i, 1, i, 4)
        i += 1
        resetbutton = QPushButton("重置")
        resetbutton.setToolTip("重置核数据为初始参数")
        grid.addWidget(resetbutton, *(i, 3))
        resetbutton.clicked.connect(self.editreset)
        checkbutton = QPushButton("确认")
        checkbutton.setToolTip("点击确认，保存修改")
        grid.addWidget(checkbutton, *(i, 4))
        checkbutton.clicked.connect(self.editok)

        self.setGeometry(230, 50, 1400, 950)
        self.show()

    def closeEvent(self, event):

        reply = QMessageBox.question(self, '消息',"退出程序？", QMessageBox.Yes | QMessageBox.No, QMessageBox.No)

        if reply == QMessageBox.Yes:
            event.accept()
        else:
            event.ignore()

    def saveData(self):
        fname = QFileDialog.getSaveFileName(self, '保存数据', '/',
                                            "Text Files (*.txt);;All Files (*)")

        if fname[0]:
            with open(fname[0], 'w') as file:
                file.write(self.Draw_Widget.dataText)

    def saveFigure(self):
        fname = QFileDialog.getSaveFileName(self, '保存图形', '/',
                                            "(*.png);;(*.eps);;(*.pdf);;(*.ps);;(*.raw);;(*.rgba);;(*.svg);;(*.svgz)")

        if fname[0]:
            plt.savefig(fname[0])

    def show_info(self):
        reply = QMessageBox.information(self, '关于', "软件名称：BurnUI\n软件类型：点燃耗计算程序\n发布时间：2020/6/15\n软件版本：1.2"
                                                    "\n新特性：\n\t1. 解决计算过程中无响应问题"
                                                    "\n\t2. 添加了进度条"
                                                    "\n-------------------------------------------------\n"
                                                    "作者：张子文\n单位：清华大学"
                                                    "工程物理系\nEmail: zhangzw17@mails.tsinghua.edu.cn"
                                                    "\n-------------------------------------------------\n"
                                                    "欢迎使用！",
                                        QMessageBox.Ok, QMessageBox.Ok)

    def help(self):
        reply = QMessageBox.information(self, '帮助', "考虑的燃耗链: U234->U235->U236 \nU238->U239->Np239->Pu239->Pu240->Pu241->Pu242\n"
                                                    "考虑的裂变产物: I135->Xe135\n"
                                                    "外步长：指整个燃耗时间段的划分，每点击一次‘计算’，程序自动计算一个外步长的时间，并停下来等待用户修改通量、功率等参数。\n"
                                                    "内步长：指燃耗方程数值求解过程中的步长，影响程序计算的精度。\n"
                                                    "定通量燃耗：指在所划分的燃耗时间段内中子通量水平保持不变。\n"
                                                    "定功率燃耗：指在所划分的燃耗时间段内反应堆功率水平保持不变，程序自动在每个内步长上根据当前核密度的变化调整中子通量。\n"
                                                    "初值：只有第一步计算时可输入各核素核密度的初值，之后每步计算都自动以上一步计算的结果作为初值。\n"
                                                    "核数据：来源于ENDF，能量0.0253eV，用户可根据需求在工具栏‘修改核数据’中进行修改。",
                                        QMessageBox.Ok, QMessageBox.Ok)

    def edit(self):
        for name, Nuc in Nuclides.items():
            Nuc.sig_a_input.setText(str(Nuc.sig_a/barn))
            Nuc.sig_g_input.setText(str(Nuc.sig_g/barn))
            Nuc.sig_f_input.setText(str(Nuc.sig_f/barn))
            Nuc.lamda_input.setText(str('%.4e' % Nuc.lamda))
            if type(Nuc) == type(I135):
                show25 = Nuc.gama25 * 100
                show49 = Nuc.gama49 * 100
                show41 = Nuc.gama41 * 100
                Nuc.gama25_input.setText(str('%.3f' % show25))
                Nuc.gama49_input.setText(str('%.3f' % show49))
                Nuc.gama41_input.setText(str('%.3f' % show41))
        self.edit_window.show()

    def editok(self):
        try:
            for name, Nuc in Nuclides.items():
                Nuc.sig_a = float(Nuc.sig_a_input.text()) * barn
                Nuc.sig_g = float(Nuc.sig_g_input.text()) * barn
                Nuc.sig_f = float(Nuc.sig_f_input.text()) * barn
                Nuc.lamda = float(Nuc.lamda_input.text())
                if type(Nuc) == type(I135):
                    Nuc.gama25 = float(Nuc.gama25_input.text()) / 100
                    Nuc.gama49 = float(Nuc.gama49_input.text()) / 100
                    Nuc.gama41 = float(Nuc.gama41_input.text()) / 100
            self.edit_window.close()
        except ValueError:
            QMessageBox.warning(self, '警告', "输入错误，请重新输入！", QMessageBox.Ok, QMessageBox.Ok)


    def editreset(self):
        reply = QMessageBox.question(self, '重置',"确认重置数据？", QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        if reply == QMessageBox.Yes:
            for name, Nuc in Nuclides.items():
                Nuc.sig_a_input.setText(str(Nuc.sig_a_origin/barn))
                Nuc.sig_g_input.setText(str(Nuc.sig_g_origin/barn))
                Nuc.sig_f_input.setText(str(Nuc.sig_f_origin/barn))
                Nuc.lamda_input.setText(str('%.4e' % Nuc.lamda_origin))
                if type(Nuc) == type(I135):
                    show25 = Nuc.gama25_origin * 100
                    show49 = Nuc.gama49_origin * 100
                    show41 = Nuc.gama41_origin * 100
                    Nuc.gama25_input.setText(str('%.3f' % show25))
                    Nuc.gama49_input.setText(str('%.3f' % show49))
                    Nuc.gama41_input.setText(str('%.3f' % show41))


# 运行程序
if __name__ == '__main__':
    app = QApplication(argv)

    # 单位初始化
    barn = 1e-24  # cm^2
    N = 1e22  # cm^-3
    minute = 60
    hr = 60 * minute
    day = 24 * hr
    year = 365 * day

    # 输入参数
    s = 0
    P = 300e6
    Dt = 3000 * day
    dt = 1 * day
    # 计算时间步
    j = int(Dt // dt)
    ti = arange(0, Dt + dt, dt)

    # 初始化核素
    Nuclides = {}
    U0 = Nuclide("U0", 0, 0)  # 空核素
    del Nuclides["U0"]  # 不需要计算空核素
    U234 = Nuclide("U234 ", 100.975, 100.908, N0=0.001, parent=U0)
    U235 = Nuclide("U235 ", 683.681, 99.3843, sig_f=586.691, N0=0.03 * 4.82, parent=U234)
    U236 = Nuclide("U236 ", 5.18107, 5.13396, parent=U235)
    U238 = Nuclide("U238 ", 2.7001772, 2.68368, N0=0.97 * 4.82, parent=U0)
    U239 = Nuclide("U239 ", 36.7557, 22.5258, T=23.5 * minute, parent=U238)
    Np239 = Nuclide("Np239", 46.0181, 0, T=2.355 * day, parent=U239)
    Pu239 = Nuclide("Pu239", 1017.53, 270.139, sig_f=747.393, parent=Np239)
    Pu240 = Nuclide("Pu240", 289.585, 289.529, parent=Pu239)
    Pu241 = Nuclide("Pu241", 1374.03, 363.047, sig_f=1012.3, T=14.4 * year, parent=Pu240)
    Pu242 = Nuclide("Pu242", 18.679, 0, parent=Pu241)
    # 裂变产物
    I135 = Product("I135 ", 0, 0, 6.386, 6.100, 7.694, parent=U0, lamda=2.87e-5)
    Xe135 = Product("Xe135", 2.7e6, 2.7e6, 0.228, 1.087, 0.255, parent=I135, lamda=2.09e-5)

    main_window = Main_window()
    app.exec()
