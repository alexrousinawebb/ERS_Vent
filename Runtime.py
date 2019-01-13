from __future__ import print_function, unicode_literals

import ODE
from Scenario import Scenario

from tqdm import tqdm
import matplotlib.pyplot as plt
import numpy as np
import os
import click
import six
from pyfiglet import figlet_format
from PyInquirer import (Token, ValidationError, Validator,
                        style_from_dict, prompt)

try:
    import colorama
    colorama.init()
except ImportError:
    colorama = None

try:
    from termcolor import colored
except ImportError:
    colored = None

style = style_from_dict({
    Token.QuestionMark: '#fac731 bold',
    Token.Answer: '#4688f1 bold',
    Token.Instruction: '',  # default
    Token.Separator: '#cc5454',
    Token.Selected: '#0abf5b',  # default
    Token.Pointer: '#673ab7 bold',
    Token.Question: '',
})

def log(string, color, font="slant", figlet=False):
    if colored:
        if not figlet:
            six.print_(colored(string, color))
        else:
            six.print_(colored(figlet_format(
                string, font=font), color))
    else:
        six.print_(string)

class Num_Validator(Validator):
    def validate(self, value):
        if len(value.text):
            try:
                float(value.text)
            except:
                raise ValidationError(
                    message="Input must be a number",
                    cursor_position=len(value.text))

            val = float(value.text)
            if val >= 0:
                return True
            else:
                raise ValidationError(
                    message="Input must be positive",
                    cursor_position=len(value.text))
        else:
            raise ValidationError(
                message="You can't leave this blank",
                cursor_position=len(value.text))

class Sensitivity():
    def sensitivity(self, scenario, value, ranges):
        for i in tqdm(ranges):
            if value == "Rupture Disc Diameter":
                scenario.D_RD = i
            elif value == "Rupture Disc Burst Pressure":
                scenario.P_RD = i
                input('press a thing')
            elif value == "Backpressure Regulator Set-Point":
                scenario.P_BPR = i
            elif value == "Hydrogen Peroxide Concentration":
                scenario.XH2O2 = i/100
            elif value == "Reactor Charge":
                scenario.mR = i
            elif value == "Contamination Factor":
                scenario.kf = i
            elif value == "Reaction Temperature":
                scenario.rxn_temp = i

            ode1 = ODE.ODE(scenario)
            ode1.initialize_heatup()
            ode1.integrate()

            if scenario.RD is True and ode1.tc == 1:
                ode1.initialize_vent(integrator='vode')
                ode1.integrate()

            self.max_P.append(ode1.max_P())
            self.max_T.append(ode1.max_T())
            self.max_conversion.append(ode1.max_conversion())
            self.max_vent.append(ode1.max_conversion())

    def plot_sensitivity(self, value, ranges):
        plt.figure(1, figsize=(10, 10))

        plt.subplot(2, 2, 1)
        plt.plot(ranges, self.max_P, color='r')
        plt.xlabel(str(value))
        plt.ylabel('Pressure (kPa)')
        plt.title("Maximum Reactor Pressure")

        plt.subplot(2, 2, 2)
        plt.plot(ranges, self.max_T, color='r')
        plt.xlabel(str(value))
        plt.ylabel('Temperature (deg C)')
        plt.title("Maximum Reactor Temperature")

        plt.subplot(2, 2, 3)
        plt.plot(ranges, self.max_conversion, color='r')
        plt.xlabel(str(value))
        plt.ylabel('Conversion (%)')
        plt.title("Maximum Reactor Conversion")

        plt.subplot(2, 2, 4)
        plt.plot(ranges, self.max_vent, color='r')
        plt.xlabel(str(value))
        plt.ylabel('Flow Rate (mol/s)')
        plt.title("Maximum Vent Flow")

        plt.show()

class Questions():
    def greeting(self):
        log("ERS Vent", color="blue", figlet=True)
        log("Welcome to XRCC Emergency Relief System Vent (Version 0.1 Beta)", "blue")

    def information(self):
        log("XRCC ERS Vent - Batch Process Emergency Relief System Modelling and Sizing Using DIERS Technology", 'blue')
        print("ERS Vent is software designed for modelling, simulation, and optimization of batch reactor "
               "emergency relief systems for reactive systems")
        print("Software is for use by registered XRCC employees only. If you believe you have recieved this "
              "software in error please contact your local administrator")
        print("This is beta software provided for use without any warranty or technical support")
        input("Press [Enter] to return...")

    def not_implemented(self):
        input("Not Yet Implemeted, Press [Enter] to return...")

    def root_menu_q(self):
        questions = [
            {
                'type': 'list',
                'name': 'Main Menu',
                'message': 'Main Menu',
                'choices': ['New Scenario', 'Load Scenario', 'Information', 'Exit'],
            },
        ]
        answers = prompt(questions, style=style)
        return answers

    def new_scenario_menu_q(self):
        questions = [
            {
                'type': 'list',
                'name': 'New Scenario',
                'message': 'New Scenario',
                'choices': [
                            'Model Scenario', 'RD/PRV Sizing', 'Sensitivity Analysis', 'Return'
                ],
            },
        ]
        answers = prompt(questions, style=style)
        return answers

    def model_scenario_menu_q(self):
        questions = [
            {
                'type': 'list',
                'name': 'Scenario',
                'message': 'Scenario Menu:',
                'choices': ['Configure Scenario', 'View Scenario Settings', 'Run Scenario', 'Return'],
            },
        ]
        answers = prompt(questions, style=style)
        return answers

    def config_scenario_menu_q(self):
        questions = [
            {
                'type': 'list',
                'name': 'Scenario Setup',
                'message': 'Scenario Configuration:',
                'choices': ['ERS Settings', 'Vessel Settings', 'Reaction Settings',
                            'PID Controller Settings (Optional)', 'Return'],
            },
        ]
        answers = prompt(questions, style=style)
        return answers

    def setup_ers_menu_design_q(self):
        questions = [
            {
                'type': 'checkbox',
                'name': 'ERS Setup',
                'message': 'Choose Model ERS Setup:',
                'choices': [
                    {
                        'name' : 'Rupture Disk (RD)'
                    },
                    {
                        'name' : 'Backpressure Regulator (BPR)'
                    }
                ],
            },
        ]
        answers = prompt(questions, style=style)
        return answers

    def setup_ers_menu_rd_q(self):
        questions = [
            {
                'type': 'input',
                'name': 'D_RD',
                'message': 'Rupture Disk Diameter (inches):',
                'validate': Num_Validator,
            },
            {
                'type': 'input',
                'name': 'P_RD',
                'message': 'Rupture Disk Burst Pressure (kPa):',
                'validate': Num_Validator,
            },
        ]
        answers = prompt(questions, style=style)
        return answers

    def setup_ers_menu_bpr_q(self):
        questions = [
            {
                'type': 'input',
                'name': 'D_BPR',
                'message': 'Backpressure Regulator Diameter (inches):',
                'default': '0.5',
                'validate': Num_Validator,
            },
            {
                'type': 'input',
                'name': 'Cv_BPR',
                'message': 'Backpressure Regulator Maximum Flow Coefficient (Cv):',
                'default': '5.5',
                'validate': Num_Validator,
            },
            {
                'type': 'input',
                'name': 'P_BPR',
                'message': 'Backpressure Regulator Setpoint (kPa):',
                'validate': Num_Validator,
            },
        ]
        answers = prompt(questions, style=style)

        return answers

    def setup_ers_menu_tf_q(self):
        questions = [
            {
                'type': 'list',
                'name': 'Two Phase?',
                'message': 'Select Emergency Relief Model:',
                'choices': ['All Vapour Venting', 'Two-Phase Bubbly', 'Two-Phase Churn-Turbulent', 'More Information'],
            },
        ]
        answers = prompt(questions, style=style)
        return answers

    def setup_ers_menu_tfinfo_q(self):
        questions = [
            {
                'type': 'list',
                'name': 'Two Phase?',
                'message': 'Select an Item to Learn More:',
                'choices': ['All Vapour Venting', 'Two-Phase Bubbly', 'Two-Phase Churn-Turbulent', 'Return'],
            },
        ]
        answers = prompt(questions, style=style)
        return answers

    def setup_vessel_menu_q(self):
        questions = [
            {
                'type': 'input',
                'name': 'VR',
                'message': 'Reactor Total Volume (gallons):',
                'validate': Num_Validator,
            },
            {
                'type': 'input',
                'name': 'AR',
                'message': 'Reactor Aspect Ratio (Height/Diameter):',
                'default': '1.5',
                'validate': Num_Validator,
            },
            {
                'type': 'input',
                'name': 'Ux',
                'message': 'Reactor Heat Transfer Coefficient (W/(m**2 K)):',
                'default': '450',
                'validate': Num_Validator,
            },
            {
                'type': 'input',
                'name': 'MAWP',
                'message': 'Vessel Maximum Allowable Working Pressure (kPa):',
                'default': '10000',
                'validate': Num_Validator,
            },
        ]
        answers = prompt(questions, style=style)
        return answers

    def setup_rxn_menu_q(self):
        questions = [
            {
                'type': 'input',
                'name': 'XH2O2',
                'message': 'Starting Hydrogen Peroxide Percentage (% w/w):',
                'default': '30',
                'validate': Num_Validator,
            },
            {
                'type': 'input',
                'name': 'mR',
                'message': 'Total Reactor Charge (kg):',
                'default': '304',
                'validate': Num_Validator,
            },
            {
                'type': 'input',
                'name': 'T0',
                'message': 'Starting Temperature (deg C):',
                'default': '25',
                'validate': Num_Validator,
            },
            {
                'type': 'input',
                'name': 'Trxn',
                'message': 'Reaction Temperature (deg C):',
                'default': '110',
                'validate': Num_Validator,
            },
            {
                'type': 'input',
                'name': 'P0',
                'message': 'Starting Headspace Pressure (kPa):',
                'default': '101.325',
                'validate': Num_Validator,
            },
            {
                'type': 'input',
                'name': 'kf',
                'message': 'Contamination Factor (1 - 10,000):',
                'default': '1',
                'validate': Num_Validator,
            },
            {
                'type': 'input',
                'name': 't_rxn',
                'message': 'Total Reaction Time (h):',
                'default': '6',
                'validate': Num_Validator,
            },
            {
                'type': 'input',
                'name': 't_cool',
                'message': 'Cooldown Time Following Reaction (h):',
                'default': '2',
                'validate': Num_Validator,
            },
        ]
        answers = prompt(questions, style=style)
        return answers

    def setup_pid_menu_q(self):
        questions = [
            {
                'type': 'input',
                'name': 'max_rate',
                'message': 'Maximum Rate of Jacket Temperature Change (deg C / min):',
                'default': '2',
                'validate': Num_Validator,
            },
            {
                'type': 'input',
                'name': 'Kp',
                'message': 'Proportional Gain (Kp):',
                'default': '0.016',
                'validate': Num_Validator,
            },
            {
                'type': 'input',
                'name': 'Ki',
                'message': 'Integral Gain (Ki):',
                'default': '0',
                'validate': Num_Validator,
            },
            {
                'type': 'input',
                'name': 'Kd',
                'message': 'Derivative Gain (Kd):',
                'default': '0',
                'validate': Num_Validator,
            },
        ]
        answers = prompt(questions, style=style)
        return answers

    def setup_plot_rt_q(self):
        questions = [
            {
                'type': 'confirm',
                'name': 'plot_rt',
                'message': 'Plot Solution in Real Time? (This will slow down the integrator)',
                'default': False,
                'validate': Num_Validator,
            }
        ]
        answers = prompt(questions, style=style)
        return answers

    def new_scenario_sensitivity_menu_q(self):
        questions = [
            {
                'type': 'list',
                'name': 'Sensitivity',
                'message': 'Sensitivity Analysis Menu',
                'choices': ['Configure Sensitivity', 'Configure Scenario', 'View Sensitivity Settings',
                            'View Scenario Settings', 'Run Sensitivity', 'Return'],
            },
        ]
        answers = prompt(questions, style=style)
        return answers

    def setup_scenario_sensitivity_config_q(self):
        questions = [
            {
                'type': 'input',
                'name': 'min',
                'message': 'Minimum Value:',
                'validate': Num_Validator,
            },
            {
                'type': 'input',
                'name': 'max',
                'message': 'Maximum Value:',
                'validate': Num_Validator,
            },
            {
                'type': 'input',
                'name': 'range',
                'message': 'Number of Data Points for Analysis:',
                'validate': Num_Validator,
            },
        ]
        answers = prompt(questions, style=style)
        return answers

    def setup_scenario_sensitivity_menu_q(self):
        questions = [
            {
                'type': 'list',
                'name': 'Sensitivity',
                'message': 'Select Factor for Sensitivity Analysis:',
                'choices': [
                    'Rupture Disc Diameter', 'Rupture Disc Burst Pressure', 'Backpressure Regulator Set-Point',
                    'Hydrogen Peroxide Concentration', 'Reactor Charge', 'Contamination Factor', 'Reaction Temperature'
                ],
            },
        ]
        answers = prompt(questions, style=style)
        return answers

class Logic(Questions, Sensitivity):

    def root_menu(self):
        while True:
            os.system('cls')
            self.greeting()

            answers = self.root_menu_q()
            if answers.get("Main Menu") == "New Scenario":
                self.new_scenario_menu()
            elif answers.get("Main Menu") == "Load Scenario":
                self.not_implemented()
            elif answers.get("Main Menu") == "Information":
                self.information()
            elif answers.get("Main Menu") == "Exit":
                print('Exiting ...')
                quit()

    def new_scenario_menu(self):
        while True:
            os.system('cls')
            self.greeting()

            answers = self.new_scenario_menu_q()
            if answers.get("New Scenario") == "Model Scenario":
                self.model_scenario_menu()
            elif answers.get("New Scenario") == "RD/PRV Sizing":
                self.not_implemented()
            elif answers.get("New Scenario") == "Sensitivity Analysis":
                self.new_scenario_sensitivity_menu()
            elif answers.get("New Scenario") == "Return":
                break

    def model_scenario_menu(self):
        while True:
            os.system('cls')
            self.greeting()

            answers = self.model_scenario_menu_q()
            if answers.get("Scenario") == "Configure Scenario":
                self.config_scenario_menu()
            elif answers.get("Scenario") == "View Scenario Settings":
                self.view_scenario_menu()
            elif answers.get("Scenario") == "Run Scenario":
                self.run_scenario_menu()
            elif answers.get("Scenario") == "Return":
                break

    def config_scenario_menu(self):
        while True:
            os.system('cls')
            self.greeting()

            answers = self.config_scenario_menu_q()
            if answers.get("Scenario Setup") == "ERS Settings":
                self.setup_ers_menu_design()
            elif answers.get("Scenario Setup") == "Vessel Settings":
                self.setup_vessel_menu()
            elif answers.get("Scenario Setup") == "Reaction Settings":
                self.setup_rxn_menu()
            elif answers.get("Scenario Setup") == "PID Controller Settings (Optional)":
                self.setup_pid_menu()
            elif answers.get("Scenario Setup") == "Return":
                break

    def setup_ers_menu_design(self):
        os.system('cls')
        self.greeting()

        answers = self.setup_ers_menu_design_q()
        if "Rupture Disk (RD)" in answers.get("ERS Setup"):
            self.RD = True
            self.TF = True
        else:
            self.RD = False

        if "Backpressure Regulator (BPR)" in answers.get("ERS Setup"):
            self.BPR = True
            self.TF = True
        else:
            self.BPR = False

        if self.RD is True:
            self.setup_ers_menu_tf()
            self.setup_ers_menu_rd()

        if self.BPR is True:
            self.setup_ers_menu_bpr()

        self.ers_stats()

    def setup_ers_menu_tf(self):
        while True:
            os.system('cls')
            self.greeting()

            answers = self.setup_ers_menu_tf_q()
            if answers.get("Two Phase?") == "All Vapour Venting":
                self.TF = False
                break
            elif answers.get("Two Phase?") == "Two-Phase Bubbly":
                self.TF = True
                self.flow_regime = 'bubbly'
                break
            elif answers.get("Two Phase?") == "Two-Phase Churn-Turbulent":
                self.TF = True
                self.flow_regime = 'churn-turbulent'
                break
            elif answers.get("Two Phase?") == "More Information":
                self.setup_ers_menu_tfinfo()

    def setup_ers_menu_tfinfo(self):
        while True:
            os.system('cls')
            self.greeting()

            answers = self.setup_ers_menu_tfinfo_q()
            if answers.get("Two Phase?") == "All Vapour Venting":
                print('Simulate reactor venting as single phase all vapor venting')
                input("Press [Enter] to continue...")
            elif answers.get("Two Phase?") == "Two-Phase Bubbly":
                print(
                    'Switch dynamically between single and two-phase ERS venting '
                    'based on the bubbly flow model.')
                print(
                    'The bubbly flow model is most appropriate for systems exhibiting foamy '
                    'or frothing behaviour.')
                input("Press [Enter] to continue...")
            elif answers.get("Two Phase?") == "Two-Phase Churn-Turbulent":
                print(
                    'Switch dynamically between single and two-phase ERS venting '
                    'based on the churn-turbulent flow model.')
                print(
                    'The churn-turbulent flow model is most appropriate for systems that do NOT '
                    'exhibit foamy or frothing behaviour.')
                input("Press [Enter] to continue...")
            elif answers.get("Two Phase?") == "Return":
                break

    def setup_ers_menu_rd(self):
        answers = self.setup_ers_menu_rd_q()
        self.D_RD = float(answers.get("D_RD"))
        self.P_RD = float(answers.get("P_RD"))

    def setup_ers_menu_bpr(self):
        answers = self.setup_ers_menu_bpr_q()
        self.D_BPR = float(answers.get("D_BPR"))
        self.P_BPR = float(answers.get("P_BPR"))
        self.BPR_max_Cv = float(answers.get("Cv_BPR"))

    def ers_stats(self):
        print(' ')
        log('ERS Settings:', 'blue')
        print('Rupture Disc:  ' + str(self.RD))
        print('Backpressure Regulator:  ' + str(self.BPR))
        if self.TF is True:
            print('Two Phase Flow:  ' + str(self.TF))
            print('Flow Regime:  ' + str(self.flow_regime))
        print(' ')

        if self.RD is True:
            log('Rupture Disc Parameters:', 'blue')
            print('Rupture Disc Diameter:  ' + str(self.D_RD) + ' in')
            print('Rupture Disc Burst Pressure:  ' + str(self.P_RD) + ' kPa')
            print(' ')

        if self.BPR is True:
            log('Backpressure Regulator Parameters:', 'blue')
            print('Backpressure Regulator Orifice Diameter:  ' + str(self.D_BPR) + ' in')
            print('Backpressure Regulator Maximum Flow Coefficient (Cv):  ' + str(self.BPR_max_Cv))
            print('Backpressure Regulator Set Point:  ' + str(self.P_BPR) + ' kPa')
            print(' ')

        input("Press [Enter] to continue...")

    def setup_vessel_menu(self):
        os.system('cls')
        self.greeting()

        answers = self.setup_vessel_menu_q()
        self.VR = float(answers.get("VR"))
        self.AR = float(answers.get("AR"))
        self.Ux = float(answers.get("Ux"))
        self.MAWP = float(answers.get("MAWP"))

        self.vessel_stats()

    def vessel_stats(self):
        print(' ')
        log('Vessel Parameters:', 'blue')
        print('Reactor Volume:  ' + str(self.VR) + ' gal')
        print('Reactor Aspect Ratio:  ' + str(self.AR))
        print('Heat Transfer Coefficient:  ' + str(self.Ux) + ' W/(m**2 K)')
        print('Maximum Allowable Working Pressure:  ' + str(self.MAWP) + ' kPa')
        print(' ')
        input("Press [Enter] to continue...")

    def setup_rxn_menu(self):
        os.system('cls')
        self.greeting()

        answers = self.setup_rxn_menu_q()
        self.XH2O2 = float(answers.get("XH2O2"))/100
        self.mR = float(answers.get("mR"))
        self.T0 = float(answers.get("T0"))
        self.rxn_temp = float(answers.get("Trxn"))
        self.P0 = float(answers.get("P0"))
        self.kf = float(answers.get("kf"))
        self.rxn_time = float(answers.get("t_rxn"))
        self.cool_time = float(answers.get("t_cool"))

        self.rxn_stats()

    def rxn_stats(self):
        print(' ')
        log('Reaction Parameters:', 'blue')
        print('Staring Hydrogen Peroxide Concentration:  ' + str(self.XH2O2*100) + ' % w/w')
        print('Starting Reactor Charge:  ' + str(self.mR) + ' kg')
        print('Starting Temperature:  ' + str(self.T0) + ' deg C')
        print('Reaction Temperature:  ' + str(self.rxn_temp) + ' deg C')
        print('Starting Headspace Pressure:  ' + str(self.P0) + ' kPa')
        print('Hydrogen Peroxide Contamination Factor:  ' + str(self.kf))
        print('Reaction Time:  ' + str(self.rxn_time) + ' h')
        print('Cooldown Time:  ' + str(self.cool_time) + ' h')
        print(' ')
        input("Press [Enter] to continue...")

    def setup_pid_menu(self):
        os.system('cls')
        self.greeting()

        answers = self.setup_pid_menu_q()
        self.max_rate = float(answers.get("max_rate"))
        self.Kp = float(answers.get("Kp"))
        self.Ki = float(answers.get("Ki"))
        self.Kd = float(answers.get("Kd"))

        self.pid_stats()

    def pid_stats(self):
        print(' ')
        log('PID Controller Configuration:', 'blue')
        print('Maximum Rate of Temperature Change in Jacket:  ' + str(self.max_rate) + ' deg C / min')
        print('Proportional Gain (Kp):  ' + str(self.Kp))
        print('Integral Gain (Ki):  ' + str(self.Ki))
        print('Derivative Gain (Kd):  ' + str(self.Kd))
        print(' ')
        input("Press [Enter] to continue...")

    def view_scenario_menu(self):

        if None in [self.VR, self.RD, self.kf]:
            print(' ')
            print('Scenario not fully specified. Update scenario configuration and try again')
            print(' ')
            input("Press [Enter] to continue...")

        else:
            self.ers_stats()
            self.vessel_stats()
            self.rxn_stats()
            self.pid_stats()

    def run_scenario_menu(self):
        if None in [self.VR, self.RD, self.kf]:
            print(' ')
            print('Scenario not fully specified. Update scenario configuration and try again')
            print(' ')
            input("Press [Enter] to continue...")

        else:
            answers = self.setup_plot_rt_q()
            plot_rt = answers.get("plot_rt")

            scen = Scenario(
                 self.VR, self.rxn_temp, self.rxn_time, self.XH2O2, self.mR, self.D_RD,
                 self.P_RD, self.P_BPR, self.D_BPR, self.BPR_max_Cv, self.P0, self.kf, self.T0,
                 self.Ux, self.AR, self.MAWP, self.max_rate,
                 self.Kp, self.Ki, self.Kd, self.flow_regime, self.BPR, self.RD, self.cool_time,
                 self.TF
                 )

            ode1 = ODE.ODE(scen)
            ode1.initialize_heatup()

            print('Integrating Heatup...')
            print(' ')

            ode1.integrate(plot_rt)
            tc = ode1.termination_code()

            print(' ')
            print(tc)
            print(' ')

            if self.RD is True and ode1.tc == 1:
                print('Integrating ERS...')
                print(' ')

                ode1.initialize_vent(integrator='vode')
                ode1.integrate(plot_rt)
                tc = ode1.termination_code()

                print(' ')
                print(tc)
                print(' ')

            self.summary_stats_menu(ode1)

    def summary_stats_menu(self, ode):
        max_P = ode.max_P()
        max_T = ode.max_T()
        max_conversion = ode.max_conversion()

        print(' ')
        log('Run Statistics:', 'blue')
        print('Maximum Pressure:  ' + str(round(max_P, 2)) + ' kPa')
        print('Maximum Temperature:  ' + str(round(max_T, 2)) + ' deg C')
        print('Maximum Conversion:  ' + str(round(max_conversion, 2)) + ' %')

        if (self.RD is True) or (self.BPR is True):
            max_vent = ode.max_vent()
            print('Maximum Vent Flowrate:  ' + str(round(max_vent, 4)) + ' mol/s')

            if self.TF is True:
                min_quality = ode.min_quality()
                print('Minimum Vent Quality:  ' + str(round(min_quality, 4)))

        ode.plot_vals()
        input("Press [Enter] to continue...")

    def new_scenario_sensitivity_menu(self):
        while True:
            os.system('cls')
            self.greeting()

            answers = self.new_scenario_sensitivity_menu_q()
            if answers.get("Sensitivity") == "Configure Sensitivity":
                self.setup_scenario_sensitivity_menu()
            elif answers.get("Sensitivity") == "Configure Scenario":
                self.config_scenario_menu()
            elif answers.get("Sensitivity") == "View Sensitivity Settings":
                self.view_sensitivity_menu()
            elif answers.get("Sensitivity") == "View Scenario Settings":
                self.view_scenario_menu()
            elif answers.get("Sensitivity") == "Run Sensitivity":
                self.run_sensitivity_menu()
            elif answers.get("Sensitivity") == "Return":
                break

    def setup_scenario_sensitivity_menu(self):
        os.system('cls')
        self.greeting()

        answers = self.setup_scenario_sensitivity_menu_q()
        self.value = answers.get("Sensitivity")
        self.setup_scenario_sensitivity_config()

    def setup_scenario_sensitivity_config(self):
        os.system('cls')
        self.greeting()

        print(' ')
        log(str(self.value), "blue")
        answers = self.setup_scenario_sensitivity_config_q()
        min = float(answers.get("min"))
        max = float(answers.get("max"))
        span = int(answers.get("range"))

        self.ranges = np.linspace(min, max, span)

    def view_sensitivity_menu(self):
        if None in [self.value]:
            print(' ')
            print('Sensitivity not specified. Update sensitivity configuration and try again')
            print(' ')
            input("Press [Enter] to continue...")

        else:
            print(' ')
            log(str(self.value) + " Sensitivity Chosen", "blue")
            print('Minimum Value:  ' + str(self.ranges[0]))
            print('Maximum Value:  ' + str(self.ranges[-1]))
            print(' ')
            input("Press [Enter] to continue...")

    def run_sensitivity_menu(self):
        if None in [self.VR, self.RD, self.kf]:
            print(' ')
            print('Scenario not fully specified. Update scenario configuration and try again')
            print(' ')
            input("Press [Enter] to continue...")

        elif None in [self.value]:
            print(' ')
            print('Sensitivity not specified. Update sensitivity configuration and try again')
            print(' ')
            input("Press [Enter] to continue...")

        else:
            scen = Scenario(
                 self.VR, self.rxn_temp, self.rxn_time, self.XH2O2, self.mR, self.D_RD,
                 self.P_RD, self.P_BPR, self.D_BPR, self.BPR_max_Cv, self.P0, self.kf, self.T0,
                 self.Ux, self.AR, self.MAWP, self.max_rate,
                 self.Kp, self.Ki, self.Kd, self.flow_regime, self.BPR, self.RD, self.cool_time,
                 self.TF
                 )

            print(' ')
            print('Starting Sensitivity Analysis')
            try:
                self.sensitivity(scen, self.value, self.ranges)
                self.plot_sensitivity(self.value, self.ranges)
            except:
                print('Something went wrong, please try again...')

class Menu(Logic):
    def __init__(self):
        #  Reactor Parameters
        self.VR = None
        self.Ux = None
        self.AR = None
        self.MAWP = None

        #  Rupture Disc Parameters
        self.D_RD = None
        self.P_RD = None

        # self.RD_params = [self.D_RD, self.P_RD]

        #  BPR Parameters
        self.D_BPR = None
        self.BPR_max_Cv = None
        self.P_BPR = None

        # self.BPR_params = [self.D_BPR, self.BPR_max_Cv, self.P_BPR]

        #  ERS Setup
        self.TF = False
        self.RD = None
        self.BPR = None
        self.flow_regime = None

        # self.ERS_params = [self.TF_vent, self.RD, self.BPR, self.flow_regime]

        #  Chemical Properties
        self.XH2O2 = None
        self.mR = None
        self.kf = None

        #  Reaction Conditions
        self.rxn_temp = None
        self.rxn_time = None
        self.T0 = None
        self.P0 = None
        self.cool_time = None

        #  PID Controller
        self.max_rate = 2
        self.Kp = 0.016
        self.Kd = 0
        self.Ki = 0

        self.max_P = []
        self.max_T = []
        self.max_conversion = []
        self.max_vent = []

        #  Sensitivity Analysis
        self.value = None
        self.ranges = None

        #  Run main program
        self.root_menu()

@click.command()
def main():
    """
    Simple CLI for sending emails using SendGrid
    """

    Menu()

if __name__ == '__main__':
    main()

# T_rxn = 110
# VR = 100
# XH2O2 = 0.3
# mR = 304
# t_rxn = 6
# D_RD = 6
# P_RD = 1000
# P_BPR = 200
#
# scen = Scenario(VR, T_rxn, t_rxn, XH2O2, mR, D_RD, P_RD, P_BPR, cooldown_time=2, kf=3000, RD=True, BPR=True,
#                 TF_vent=True)
#
# # optimize = ODE.Prog(scen)
# # D_RD = optimize.vent_opt()
# # print(D_RD)
#
# ode1 = ODE.ODE(scen)
# ode1.initialize_heatup()
# ode1.integrate(plot_rt=False)
# ode1.initialize_vent(integrator='vode')
# ode1.integrate(plot_rt=False)
# print(ode1.max_P())
# ode1.plot_vals()
