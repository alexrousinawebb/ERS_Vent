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

def menu_a():
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

def menu_newscen():
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

def menu_model_a():
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

def main_menu():
    os.system('cls')

    while True:

        log("ERS Vent", color="blue", figlet=True)
        log("Welcome to XRCC Emergency Relief System Vent Version 0.1 Beta", "blue")

        main_m = menu_a()

        if main_m.get("Main Menu") == "New Scenario":
            new_scenario()
            os.system('cls')
        elif main_m.get("Main Menu") == "Load Scenario":
            pprint('i will load a scenario')
        elif main_m.get("Main Menu") == "Information":
            information()
            os.system('cls')
        elif main_m.get("Main Menu") == "Exit":
            print('Exiting ...')
            break

def new_scenario():
    os.system('cls')

    while True:
        log("ERS Vent", color="blue", figlet=True)
        log("Welcome to XRCC Emergency Relief System Vent Version 0.1 Beta", "blue")

        newscen_m = menu_newscen()

        if newscen_m.get("New Scenario") == "Model Scenario":
            new_scenario()
            os.system('cls')
        elif newscen_m.get("New Scenario") == "RD/PRV Sizing":
            input("Not Yet Implemeted, Press [Enter] to return...")
            os.system('cls')
        elif newscen_m.get("New Scenario") == "Sensitivity Analysis":
            input("Not Yet Implemeted, Press [Enter] to return...")
            os.system('cls')
        elif newscen_m.get("New Scenario") == "Return":
            break

def model_scenario():
    os.system('cls')

    while True:
        newscen_m = menu_model_a()

        if newscen_m.get("Venting") is True:
            RD = True
            new_scenario()
            os.system('cls')
        else:
            RD = False


def information():
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

    main_menu()

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

