from __future__ import print_function, unicode_literals

import ODE
from Scenario import Scenario
import os
import click
import six
from pyfiglet import figlet_format
from pprint import pprint
from PyInquirer import (Token, ValidationError, Validator, print_json, prompt,
                        style_from_dict, prompt, Separator)

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

class Questions():
    def menu_a(self):
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

    def menu_newscen(self):
        questions = [
            {
                'type': 'list',
                'name': 'New Scenario',
                'message': 'New Scenario',
                'choices': ['Model Scenario', 'RD/PRV Sizing', 'Sensitivity Analysis', 'Return'],
            },
        ]
        answers = prompt(questions, style=style)
        return answers

    def menu_model_a(self):
        questions = [
            {
                'type': 'confirm',
                'name': 'Venting?',
                'message': 'Simulate reactor emergency relief venting?',
                'default': True
            },
        ]
        answers = prompt(questions, style=style)
        return answers

    def menu_model_b(self):
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

    def menu_model_b_info(self):
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

    def menu_setup_a(self):
        questions = [
            {
                'type': 'list',
                'name': 'Scenario',
                'message': 'Scenario Menu:',
                'choices': ['Configure Scenario', 'View Scenario Settings', 'Run Scenario', 'Exit'],
            },
        ]
        answers = prompt(questions, style=style)
        return answers

    def menu_setup_b(self):
        questions = [
            {
                'type': 'list',
                'name': 'Scenario Setup',
                'message': 'Scenario Configuration:',
                'choices': ['ERS Settings', 'Vessel Settings', 'Reaction Settings', 'Return'],
            },
        ]
        answers = prompt(questions, style=style)
        return answers

    def menu_setup_c(self):
        questions = [
            {
                'type': 'list',
                'name': 'ERS Setup',
                'message': 'Configure Model ERS Settings:',
                'choices': ['ERS Design', 'Vessel Settings', 'Reaction Settings', 'Return'],
            },
        ]
        answers = prompt(questions, style=style)

        return answers
    def menu_setup_d(self):
        questions = [
            {
                'type': 'list',
                'name': 'Configure Reaction Vessel Settings',
                'message': 'Scenario Configuration:',
                'choices': ['ERS Settings', 'Vessel Settings', 'Reaction Settings', 'Return'],
            },
        ]
        answers = prompt(questions, style=style)
        return answers

    def menu_setup_e(self):
        questions = [
            {
                'type': 'list',
                'name': 'Configure Reaction Kinetics',
                'message': 'Scenario Configuration:',
                'choices': ['ERS Settings', 'Vessel Settings', 'Reaction Settings', 'Return'],
            },
        ]
        answers = prompt(questions, style=style)
        return answers

class Menu(Questions):
    def __init__(self):
        #  Reactor Parameters
        self.VR = None
        self.Ux = None
        self.AR = None
        self.MAWP = None

        #  ERS Design Parameters
        self.RD = None
        self.BPR = None
        self.D_RD = None
        self.D_BPR = None
        self.BPR_max_Cv = None
        self.P_RD = None
        self.P_BPR = None
        self.flow_regime = None
        self.TF_vent = None

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
        self.max_rate = None
        self.Kp = None
        self.Kd = None
        self.Ki = None

        self.main_menu()

    def greeting(self):
        log("ERS Vent", color="blue", figlet=True)
        log("Welcome to XRCC Emergency Relief System Vent (Version 0.1 Beta)", "blue")

    def main_menu(self):
        os.system('cls')

        while True:

            self.greeting()

            main_m = self.menu_a()

            if main_m.get("Main Menu") == "New Scenario":
                self.new_scenario()
                os.system('cls')
            elif main_m.get("Main Menu") == "Load Scenario":
                input("Not Yet Implemeted, Press [Enter] to return...")
                os.system('cls')
            elif main_m.get("Main Menu") == "Information":
                self.information()
                os.system('cls')
            elif main_m.get("Main Menu") == "Exit":
                print('Exiting ...')
                quit()

    def new_scenario(self):
        os.system('cls')

        while True:
            self.greeting()

            newscen_m = self.menu_newscen()

            if newscen_m.get("New Scenario") == "Model Scenario":
                self.model_scenario()
                os.system('cls')
            elif newscen_m.get("New Scenario") == "RD/PRV Sizing":
                input("Not Yet Implemeted, Press [Enter] to return...")
                os.system('cls')
            elif newscen_m.get("New Scenario") == "Sensitivity Analysis":
                input("Not Yet Implemeted, Press [Enter] to return...")
                os.system('cls')
            elif newscen_m.get("New Scenario") == "Return":
                break

    def model_scenario(self):
        os.system('cls')

        while True:
            self.greeting()

            newscen_m = self.menu_model_a()

            if newscen_m.get("Venting?") is True:
                self.RD = True

                while True:
                    os.system('cls')
                    self.greeting()

                    newscen_v = self.menu_model_b()

                    if newscen_v.get("Two Phase?") == "All Vapour Venting":
                        self.TF_vent = False
                        self.setup_main()
                    elif newscen_v.get("Two Phase?") == "Two-Phase Bubbly":
                        self.TF_vent = True
                        self.flow_regime = 'bubbly'
                        self.setup_main()
                    elif newscen_v.get("Two Phase?") == "Two-Phase Churn-Turbulent":
                        self.TF_vent = True
                        self.flow_regime = 'churn-turbulent'
                        self.setup_main()
                    elif newscen_v.get("Two Phase?") == "More Information":
                        os.system('cls')

                        while True:
                            self.greeting()

                            newscen_i = self.menu_model_b_info()

                            if newscen_i.get("Two Phase?") == "All Vapour Venting":
                                print('Simulate reactor venting as single phase all vapor venting')
                                input("Press [Enter] to continue...")
                                os.system('cls')
                            elif newscen_i.get("Two Phase?") == "Two-Phase Bubbly":
                                print(
                                    'Switch dynamically between single and two-phase ERS venting '
                                    'based on the bubbly flow model.')
                                print(
                                    'The bubbly flow model is most appropriate for systems exhibiting foamy '
                                    'or frothing behaviour.')
                                input("Press [Enter] to continue...")
                                os.system('cls')
                            elif newscen_i.get("Two Phase?") == "Two-Phase Churn-Turbulent":
                                print(
                                    'Switch dynamically between single and two-phase ERS venting '
                                    'based on the churn-turbulent flow model.')
                                print(
                                    'The churn-turbulent flow model is most appropriate for systems that do NOT '
                                    'exhibit foamy or frothing behaviour.')
                                input("Press [Enter] to continue...")
                                os.system('cls')
                            elif newscen_i.get("Two Phase?") == "Return":
                                break
                else:
                    self.RD = False
                    self.setup_main()

    def setup_main(self):
        os.system('cls')

        while True:
            self.greeting()

            newscen_m = self.menu_setup_a()

            if newscen_m.get("Scenario") == "Configure Scenario":
                self.setup_scenario()
                os.system('cls')
            elif newscen_m.get("Scenario") == "View Scenario Settings":
                self.view_scenario()
                os.system('cls')
            elif newscen_m.get("Scenario") == "Run Scenario":
                self.run_scenario()
                os.system('cls')
            elif newscen_m.get("Scenario") == "Exit":
                print('Exiting ...')
                quit()

    def setup_scenario(self):
        os.system('cls')

        while True:
            self.greeting()

            newscen_m = self.menu_setup_b()

            if newscen_m.get("Scenario Setup") == "ERS Settings":
                self.setup_ERS()
                os.system('cls')
            elif newscen_m.get("Scenario Setup") == "Vessel Settings":
                self.setup_vessel()
                os.system('cls')
            elif newscen_m.get("Scenario Setup") == "Reaction Settings":
                self.setup_rxn()
                os.system('cls')
            elif newscen_m.get("Scenario Setup") == "Return":
                break

    def setup_ERS(self):
        os.system('cls')

        while True:
            self.greeting()

            newscen_m = self.menu_setup_c()

    def setup_vessel(self):
        os.system('cls')

        while True:
            self.greeting()

            newscen_m = self.menu_setup_d()

    def setup_rxn(self):
        os.system('cls')

        while True:
            self.greeting()

            newscen_m = self.menu_setup_e()

    def view_scenario(self):
        os.system('cls')

        while True:
            self.greeting()

            newscen_m = self.menu_setup_a()

            if newscen_m.get("Scenario") == "Configure Scenario":
                self.model_scenario()
                os.system('cls')
            elif newscen_m.get("New Scenario") == "View Scenario Settings":
                input("Not Yet Implemeted, Press [Enter] to return...")
                os.system('cls')
            elif newscen_m.get("New Scenario") == "Run Scenario":
                input("Not Yet Implemeted, Press [Enter] to return...")
                os.system('cls')
            elif newscen_m.get("New Scenario") == "Return":
                break

    def run_scenario(self):
        os.system('cls')

        while True:
            self.greeting()

            newscen_m = self.menu_setup_a()

            if newscen_m.get("Scenario") == "Configure Scenario":
                self.model_scenario()
                os.system('cls')
            elif newscen_m.get("New Scenario") == "View Scenario Settings":
                input("Not Yet Implemeted, Press [Enter] to return...")
                os.system('cls')
            elif newscen_m.get("New Scenario") == "Run Scenario":
                input("Not Yet Implemeted, Press [Enter] to return...")
                os.system('cls')
            elif newscen_m.get("New Scenario") == "Return":
                break

    def information(self):
        log("XRCC ERS Vent - Batch Process Emergency Relief System Modelling and Sizing Using DIERS Technology", 'blue')
        print("ERS Vent is software designed for modelling, simulation, and optimization of batch reactor "
               "emergency relief systems for reactive systems")
        print("Software is for use by registered XRCC employees only. If you believe you have recieved this "
              "software in error please contact your local administrator")
        print("This is beta software provided for use without any warranty or technical support")
        input("Press [Enter] to return...")

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

