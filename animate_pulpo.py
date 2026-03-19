#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Animated Octopus for Snakemake
"""

import sys
import time
import threading

# Full image of the octopus
octopus = """

PULPO V.1.0 is running...
⠀⠀⠀⠀⠀⠀⢀⣀⣠⣀⣀⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⣠⣾⣿⣿⣿⣿⣿⣿⣷⣦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⢠⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣷⡀⠀⠀⠀⣠⣶⣾⣷⣶⣄⠀⠀⠀⠀⠀
⠀⠀⠀⢸⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣧⠀⠀⢰⣿⠟⠉⠻⣿⣿⣷⠀⠀⠀⠀
⠀⠀⠀⠈⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⠿⢷⣄⠘⠿⠀⠀⠀⢸⣿⣿⡆⠀⠀⠀
⠀⠀⠀⠀⠈⠿⣿⣿⣿⣿⣿⣀⣸⣿⣷⣤⣴⠟⠀⠀⠀⠀⢀⣼⣿⣿⠁⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠈⠙⣛⣿⣿⣿⣿⣿⣿⣿⣿⣦⣀⣀⣀⣴⣾⣿⣿⡟⠀⠀⠀⠀
⠀⠀⠀⢀⣠⣴⣾⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡿⠟⠋⣠⣤⣀⠀⠀
⠀⠀⣴⣿⣿⣿⠿⠟⠛⠛⢛⣿⣿⣿⣿⣿⣿⣧⡈⠉⠁⠀⠀⠀⠈⠉⢻⣿⣧⠀
⠀⣼⣿⣿⠋⠀⠀⠀⠀⢠⣾⣿⣿⠟⠉⠻⣿⣿⣿⣦⣄⠀⠀⠀⠀⠀⣸⣿⣿⠃
⠀⣿⣿⡇⠀⠀⠀⠀⠀⣿⣿⡿⠃⠀⠀⠀⠈⠛⢿⣿⣿⣿⣿⣶⣿⣿⣿⡿⠋⠀
⠀⢿⣿⣧⡀⠀⣶⣄⠘⣿⣿⡇⠀⠀⠠⠶⣿⣶⡄⠈⠙⠛⠻⠟⠛⠛⠁⠀⠀⠀
⠀⠈⠻⣿⣿⣿⣿⠏⠀⢻⣿⣿⣄⠀⠀⠀⣸⣿⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠻⣿⣿⣿⣶⣾⣿⣿⠃⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠙⠛⠛⠛⠋⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀

Thank you for using us!

"""

running = True

def animate():
    for _ in range(5): 
        if not running:
            break
        sys.stdout.write("\033[H\033[J")  
        print(octopus)
        print("\n🐙 Execution in progress... PULPO is running...")
        time.sleep(0.5)

def start_animation():
    """Inicia la animación en un hilo"""
    global running
    running = True
    thread = threading.Thread(target=animate)
    thread.daemon = True
    thread.start()
    return thread

def stop_animation():
    """Detiene la animación y muestra el resultado final"""
    global running
    running = False
    sys.stdout.write("\033[H\033[J")  # Limpiar pantalla
    if len(sys.argv) > 1 and sys.argv[1] == "success":
        print("✅🎉 PULPO successfully completed!\n")
    else:
        print("❌💥 PULPO ended with errors!\n")

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "start":
        start_animation().join()
    else:
        stop_animation()
