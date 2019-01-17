"""
RUNTIME
===========================================================================
Module containing all command line interface backend.
"""

from __future__ import print_function, unicode_literals

import os
import string

import click
import dill
import matplotlib.pyplot as plt
import numpy as np
import six
from PyInquirer import (Token, ValidationError, Validator,
                        style_from_dict, prompt)
from prettytable import PrettyTable
from pyfiglet import figlet_format
from tqdm import tqdm

import ODE
from Scenario import Scenario

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

class Name_Validator(Validator):
    def validate(self, value):
        invalidChars = set(string.punctuation.replace("_", ""))
        if len(value.text):
            if any(char in invalidChars for char in value.text):
                raise ValidationError(
                    message="Input cannot contain special characters",
                    cursor_position=len(value.text))
            else:
                if ' ' in value.text:
                    raise ValidationError(
                        message="Input cannot contain spaces",
                        cursor_position=len(value.text))
                else:
                    return True
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

            self.data.max_P.append(ode1.max_P())
            self.data.max_T.append(ode1.max_T())
            self.data.max_conversion.append(ode1.max_conversion())
            self.data.max_vent.append(ode1.max_conversion())

    def plot_sensitivity(self, value, ranges):
        plt.figure(1, figsize=(10, 10))

        plt.subplot(2, 2, 1)
        plt.plot(ranges, self.data.max_P, color='r')
        plt.xlabel(str(value))
        plt.ylabel('Pressure (kPa)')
        plt.title("Maximum Reactor Pressure")

        plt.subplot(2, 2, 2)
        plt.plot(ranges, self.data.max_T, color='r')
        plt.xlabel(str(value))
        plt.ylabel('Temperature (deg C)')
        plt.title("Maximum Reactor Temperature")

        plt.subplot(2, 2, 3)
        plt.plot(ranges, self.data.max_conversion, color='r')
        plt.xlabel(str(value))
        plt.ylabel('Conversion (%)')
        plt.title("Maximum Reactor Conversion")

        plt.subplot(2, 2, 4)
        plt.plot(ranges, self.data.max_vent, color='r')
        plt.xlabel(str(value))
        plt.ylabel('Flow Rate (mol/s)')
        plt.title("Maximum Vent Flow")

        plt.show()

    def table_sensitivity(self, value, ranges):
        t = PrettyTable()
        column_names = [str(value), 'Maximum Pressure (kPa)', 'Maximum Temperature (deg C)', 'Maximum Conversion (%)',
                         'Maximum Vent Flow Rate (mol/s)']
        t.add_column(column_names[0], [round(i, 2) for i in ranges])
        t.add_column(column_names[1], [round(i, 2) for i in self.data.max_P])
        t.add_column(column_names[2], [round(i, 2) for i in self.data.max_T])
        t.add_column(column_names[3], [round(i, 2) for i in self.data.max_conversion])
        t.add_column(column_names[4], [round(i, 2) for i in self.data.max_vent])
        print(t)

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
                'message': 'New Scenario: ' + str(self.data.name),
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
                'message': str(self.data.name) + ' Scenario Menu:',
                'choices': ['Configure Scenario', 'View Scenario Settings', 'Run Scenario', 'View Scenario Results',
                            'Save Data', 'Return'],
            },
        ]
        answers = prompt(questions, style=style)
        return answers

    def config_scenario_menu_q(self):
        questions = [
            {
                'type': 'list',
                'name': 'Scenario Setup',
                'message': str(self.data.name) +' Scenario Configuration:',
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
                'message': str(self.data.name) + ' Sensitivity Analysis Menu',
                'choices': ['Configure Sensitivity', 'Configure Scenario', 'View Sensitivity Settings',
                            'View Scenario Settings', 'Run Sensitivity', 'View Sensitivity Results', 'Save Data',
                            'Return'],
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

    def input_scenario_name_q(self):
        questions = [
            {
                'type': 'input',
                'name': 'scenario_name',
                'message': 'Input a name for this scenario:',
                'validate': Name_Validator,
            },
        ]
        answers = prompt(questions, style=style)
        return answers

    def load_data_q(self, files):
        questions = [
            {
                'type': 'list',
                'name': 'files',
                'message': 'Select data file to load:',
                'choices': files,
            },
        ]
        answers = prompt(questions, style=style)
        return answers

class Runtime(Questions, Sensitivity):
    def __init__(self):
        self.data = None

        self.root_menu()

    def root_menu(self):
        while True:
            os.system('cls')
            self.greeting()

            answers = self.root_menu_q()
            if answers.get("Main Menu") == "New Scenario":
                self.input_scenario_name()
                self.new_scenario_menu()
            elif answers.get("Main Menu") == "Load Scenario":
                self.load_data()
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
            elif answers.get("Scenario") == "View Scenario Results":
                self.summary_stats_menu()
            elif answers.get("Scenario") == "Save Data":
                self.save_data()
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
            self.data.RD = True
            self.data.TF = True
        else:
            self.data.RD = False

        if "Backpressure Regulator (BPR)" in answers.get("ERS Setup"):
            self.data.BPR = True
            self.data.TF = True
        else:
            self.data.BPR = False

        if self.data.RD is True:
            self.setup_ers_menu_tf()
            self.setup_ers_menu_rd()

        if self.data.BPR is True:
            self.setup_ers_menu_bpr()

        self.ers_stats()

    def setup_ers_menu_tf(self):
        while True:
            os.system('cls')
            self.greeting()

            answers = self.setup_ers_menu_tf_q()
            if answers.get("Two Phase?") == "All Vapour Venting":
                self.data.TF = False
                break
            elif answers.get("Two Phase?") == "Two-Phase Bubbly":
                self.data.TF = True
                self.data.flow_regime = 'bubbly'
                break
            elif answers.get("Two Phase?") == "Two-Phase Churn-Turbulent":
                self.data.TF = True
                self.data.flow_regime = 'churn-turbulent'
                break
            elif answers.get("Two Phase?") == "More Information":
                self.data.setup_ers_menu_tfinfo()

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
        self.data.D_RD = float(answers.get("D_RD"))
        self.data.P_RD = float(answers.get("P_RD"))

    def setup_ers_menu_bpr(self):
        answers = self.setup_ers_menu_bpr_q()
        self.data.D_BPR = float(answers.get("D_BPR"))
        self.data.P_BPR = float(answers.get("P_BPR"))
        self.data.BPR_max_Cv = float(answers.get("Cv_BPR"))

    def ers_stats(self):
        print(' ')
        log('ERS Settings:', 'blue')
        print('Rupture Disc:  ' + str(self.data.RD))
        print('Backpressure Regulator:  ' + str(self.data.BPR))
        if self.data.TF is True:
            print('Two Phase Flow:  ' + str(self.data.TF))
            print('Flow Regime:  ' + str(self.data.flow_regime))
        print(' ')

        if self.data.RD is True:
            log('Rupture Disc Parameters:', 'blue')
            print('Rupture Disc Diameter:  ' + str(self.data.D_RD) + ' in')
            print('Rupture Disc Burst Pressure:  ' + str(self.data.P_RD) + ' kPa')
            print(' ')

        if self.data.BPR is True:
            log('Backpressure Regulator Parameters:', 'blue')
            print('Backpressure Regulator Orifice Diameter:  ' + str(self.data.D_BPR) + ' in')
            print('Backpressure Regulator Maximum Flow Coefficient (Cv):  ' + str(self.data.BPR_max_Cv))
            print('Backpressure Regulator Set Point:  ' + str(self.data.P_BPR) + ' kPa')
            print(' ')

        input("Press [Enter] to continue...")

    def setup_vessel_menu(self):
        os.system('cls')
        self.greeting()

        answers = self.setup_vessel_menu_q()
        self.data.VR = float(answers.get("VR"))
        self.data.AR = float(answers.get("AR"))
        self.data.Ux = float(answers.get("Ux"))
        self.data.MAWP = float(answers.get("MAWP"))

        self.vessel_stats()

    def vessel_stats(self):
        print(' ')
        log('Vessel Parameters:', 'blue')
        print('Reactor Volume:  ' + str(self.data.VR) + ' gal')
        print('Reactor Aspect Ratio:  ' + str(self.data.AR))
        print('Heat Transfer Coefficient:  ' + str(self.data.Ux) + ' W/(m**2 K)')
        print('Maximum Allowable Working Pressure:  ' + str(self.data.MAWP) + ' kPa')
        print(' ')
        input("Press [Enter] to continue...")

    def setup_rxn_menu(self):
        os.system('cls')
        self.greeting()

        answers = self.setup_rxn_menu_q()
        self.data.XH2O2 = float(answers.get("XH2O2"))/100
        self.data.mR = float(answers.get("mR"))
        self.data.T0 = float(answers.get("T0"))
        self.data.rxn_temp = float(answers.get("Trxn"))
        self.data.P0 = float(answers.get("P0"))
        self.data.kf = float(answers.get("kf"))
        self.data.rxn_time = float(answers.get("t_rxn"))
        self.data.cool_time = float(answers.get("t_cool"))

        self.rxn_stats()

    def rxn_stats(self):
        print(' ')
        log('Reaction Parameters:', 'blue')
        print('Staring Hydrogen Peroxide Concentration:  ' + str(self.data.XH2O2*100) + ' % w/w')
        print('Starting Reactor Charge:  ' + str(self.data.mR) + ' kg')
        print('Starting Temperature:  ' + str(self.data.T0) + ' deg C')
        print('Reaction Temperature:  ' + str(self.data.rxn_temp) + ' deg C')
        print('Starting Headspace Pressure:  ' + str(self.data.P0) + ' kPa')
        print('Hydrogen Peroxide Contamination Factor:  ' + str(self.data.kf))
        print('Reaction Time:  ' + str(self.data.rxn_time) + ' h')
        print('Cooldown Time:  ' + str(self.data.cool_time) + ' h')
        print(' ')
        input("Press [Enter] to continue...")

    def setup_pid_menu(self):
        os.system('cls')
        self.greeting()

        answers = self.setup_pid_menu_q()
        self.data.max_rate = float(answers.get("max_rate"))
        self.data.Kp = float(answers.get("Kp"))
        self.data.Ki = float(answers.get("Ki"))
        self.data.Kd = float(answers.get("Kd"))

        self.pid_stats()

    def pid_stats(self):
        print(' ')
        log('PID Controller Configuration:', 'blue')
        print('Maximum Rate of Temperature Change in Jacket:  ' + str(self.data.max_rate) + ' deg C / min')
        print('Proportional Gain (Kp):  ' + str(self.data.Kp))
        print('Integral Gain (Ki):  ' + str(self.data.Ki))
        print('Derivative Gain (Kd):  ' + str(self.data.Kd))
        print(' ')
        input("Press [Enter] to continue...")

    def view_scenario_menu(self):

        if None in [self.data.VR, self.data.RD, self.data.kf]:
            print(' ')
            print('Scenario not fully specified. Update scenario configuration and try again')
            print(' ')
            input("Press [Enter] to continue...")

        else:
            print('Scenario stats for ' + str(self.data.name))
            self.ers_stats()
            self.vessel_stats()
            self.rxn_stats()
            self.pid_stats()

    def run_scenario_menu(self):
        if None in [self.data.VR, self.data.RD, self.data.kf]:
            print(' ')
            print('Scenario not fully specified. Update scenario configuration and try again')
            print(' ')
            input("Press [Enter] to continue...")

        else:
            answers = self.setup_plot_rt_q()
            plot_rt = answers.get("plot_rt")

            scen = Scenario(
                 self.data.VR, self.data.rxn_temp, self.data.rxn_time, self.data.XH2O2, self.data.mR, self.data.D_RD,
                 self.data.P_RD, self.data.P_BPR, self.data.D_BPR, self.data.BPR_max_Cv, self.data.P0, self.data.kf,
                 self.data.T0, self.data.Ux, self.data.AR, self.data.MAWP, self.data.max_rate,
                 self.data.Kp, self.data.Ki, self.data.Kd, self.data.flow_regime, self.data.BPR, self.data.RD,
                 self.data.cool_time, self.data.TF
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

            if self.data.RD is True and ode1.tc == 1:
                print('Integrating ERS...')
                print(' ')

                ode1.initialize_vent(integrator='vode')
                ode1.integrate(plot_rt)
                tc = ode1.termination_code()

                print(' ')
                print(tc)
                print(' ')

            self.data.ode = ode1
            self.summary_stats_menu()

    def summary_stats_menu(self):
        if self.data.ode is None:
            print(' ')
            print('No data to display.')
            print(' ')
            input("Press [Enter] to continue...")
        else:
            max_P = self.data.ode.max_P()
            max_T = self.data.ode.max_T()
            max_conversion = self.data.ode.max_conversion()

            print(' ')
            log('Run Statistics:', 'blue')
            print('Maximum Pressure:  ' + str(round(max_P, 2)) + ' kPa')
            print('Maximum Temperature:  ' + str(round(max_T, 2)) + ' deg C')
            print('Maximum Conversion:  ' + str(round(max_conversion, 2)) + ' %')

            if (self.data.RD is True) or (self.data.BPR is True):
                max_vent = self.data.ode.max_vent()
                print('Maximum Vent Flowrate:  ' + str(round(max_vent, 4)) + ' mol/s')

                if self.data.TF is True:
                    min_quality = self.data.ode.min_quality()
                    print('Minimum Vent Quality:  ' + str(round(min_quality, 4)))

            self.data.ode.plot_vals()
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
            elif answers.get("Sensitivity") == "View Sensitivity Results":
                self.view_sensitivity_results_menu()
            elif answers.get("Sensitivity") == "Save Data":
                self.save_data()
            elif answers.get("Sensitivity") == "Return":
                break

    def setup_scenario_sensitivity_menu(self):
        os.system('cls')
        self.greeting()

        answers = self.setup_scenario_sensitivity_menu_q()
        self.data.value = answers.get("Sensitivity")
        self.setup_scenario_sensitivity_config()

    def setup_scenario_sensitivity_config(self):
        os.system('cls')
        self.greeting()

        print(' ')
        log(str(self.data.value), "blue")
        answers = self.setup_scenario_sensitivity_config_q()
        min = float(answers.get("min"))
        max = float(answers.get("max"))
        span = int(answers.get("range"))

        self.data.ranges = np.linspace(min, max, span)

    def view_sensitivity_menu(self):
        if None in [self.data.value]:
            print(' ')
            print('Sensitivity not specified. Update sensitivity configuration and try again')
            print(' ')
            input("Press [Enter] to continue...")

        else:
            print(' ')
            log(str(self.data.value) + " Sensitivity Chosen", "blue")
            print('Minimum Value:  ' + str(self.data.ranges[0]))
            print('Maximum Value:  ' + str(self.data.ranges[-1]))
            print(' ')
            input("Press [Enter] to continue...")

    def view_sensitivity_results_menu(self):
        if len(self.data.max_P) == 0:
            print(' ')
            print('Sensitivity not yet calculted. Run sensitivity and try again')
            print(' ')
            input("Press [Enter] to continue...")
        else:
            self.stats_sensitivity()
            print(' ')
            input("Press [Enter] to continue...")

    def run_sensitivity_menu(self):
        if None in [self.data.VR, self.data.RD, self.data.kf]:
            print(' ')
            print('Scenario not fully specified. Update scenario configuration and try again')
            print(' ')
            input("Press [Enter] to continue...")

        elif None in [self.data.value]:
            print(' ')
            print('Sensitivity not specified. Update sensitivity configuration and try again')
            print(' ')
            input("Press [Enter] to continue...")

        else:
            scen = Scenario(
                 self.data.VR, self.data.rxn_temp, self.data.rxn_time, self.data.XH2O2, self.data.mR, self.data.D_RD,
                 self.data.P_RD, self.data.P_BPR, self.data.D_BPR, self.data.BPR_max_Cv, self.data.P0, self.data.kf,
                 self.data.T0, self.data.Ux, self.data.AR, self.data.MAWP, self.data.max_rate,
                 self.data.Kp, self.data.Ki, self.data.Kd, self.data.flow_regime, self.data.BPR, self.data.RD,
                 self.data.cool_time, self.data.TF
                 )

            print(' ')
            print('Starting Sensitivity Analysis for ' + str(self.data.name))
            try:
                self.sensitivity(scen, self.data.value, self.data.ranges)
                self.stats_sensitivity()
            except:
                print('Something went wrong, please try again...')

            print(' ')
            input("Press [Enter] to continue...")

    def stats_sensitivity(self):
        print(' ')
        log(str(self.data.name) + " Sensitivity Summary: " + str(self.data.value), "blue")
        print(' ')
        self.table_sensitivity(self.data.value, self.data.ranges)
        self.plot_sensitivity(self.data.value, self.data.ranges)

    def save_data(self):
        if None in [self.data.VR, self.data.RD, self.data.kf]:
            print(' ')
            print('Scenario data is incomplete, cannot save.')
            print(' ')
            input("Press [Enter] to continue...")
        else:
            try:
                with open(str(self.data.name)+'.vent', 'wb') as f:
                    dill.dump(self.data, f)
                print(' ')
                print('Data saved successfully.')
                print(' ')
                input("Press [Enter] to continue...")
            except:
                print(' ')
                print('Something went wrong, please try again.')
                print(' ')
                input("Press [Enter] to continue...")

    def load_data(self):
        directory = [f for f in os.listdir() if f.endswith('.vent')]
        if len(directory) == 0:
            print(' ')
            print('No data files present in root directory.')
            print(' ')
            input("Press [Enter] to continue...")
        else:
            answers = self.load_data_q(directory)
            file_name = answers.get("files")
            try:
                with open(file_name, 'rb') as f:
                    self.data = dill.load(f)
                print(' ')
                print('Session loaded successfully')
                print(' ')
                input("Press [Enter] to continue...")
                self.new_scenario_menu()
            except:
                print(' ')
                print('File could not be loaded.')
                print(' ')
                input("Press [Enter] to continue...")

    def input_scenario_name(self):
        self.data = Data()
        answers = self.input_scenario_name_q()
        self.data.name = answers.get("scenario_name")

class Data():
    def __init__(self):
        self.ode = None

        #  Reactor Parameters
        self.VR = None
        self.Ux = None
        self.AR = None
        self.MAWP = None

        #  Rupture Disc Parameters
        self.D_RD = None
        self.P_RD = None

        #  BPR Parameters
        self.D_BPR = None
        self.BPR_max_Cv = None
        self.P_BPR = None

        #  ERS Setup
        self.TF = False
        self.RD = None
        self.BPR = None
        self.flow_regime = None

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

        #  Sensitivity Analysis
        self.value = None
        self.ranges = None
        self.max_P = []
        self.max_T = []
        self.max_conversion = []
        self.max_vent = []

        #  Misc
        self.name = None

        #  Run main program


@click.command()
def main():
    """
    CLI for ERS Vent
    """

    Runtime()

if __name__ == '__main__':
    main()
