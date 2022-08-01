import tkinter as tk
from tkinter import *
from tkinter.ttk import *
import time
from turtle import update
import numpy as np
from PIL import ImageTk, Image
import PIL.Image
import os
import threading

"""def some_function():
    print('Do stuff')
    time.sleep(3)

if __name__ == '__main__':
    # Create window
    root = tk.Tk()
    #root.geometry('200x150')
    root.eval('tk::PlaceWindow . center')
    root.configure(bg='white')
    conv = True
    label_var = tk.StringVar()
    label = tk.Label(root, textvariable=label_var, font=['Arial',50], fg='#0059b3')
    
    label_var.set('\n -------------------------- \n  DIDO \n -------------------------- \n ')
    label.pack()
    def run():
        label_var.set(np.random.rand(1))
        return root.after(300,run)
        # Create label
    progress = Progressbar(root, orient = HORIZONTAL, length = 100, mode = 'determinate')
    label.pack()

            #a,b = runpat()
    

    root.after(800,run)
    root.mainloop()"""
    
    # Update and show window once
   
class clinicwindow:
    def __init__(self, root):
        self.root = root
        self.conv = None
    
    def setwin(self):
        self.root.geometry('660x440')

        PictureLabel= tk.Label(self.root)
        PictureLabel.pack()
        self.img = ImageTk.PhotoImage(file='code_development\\Photonics.jpg')
        PictureLabel.configure(image=self.img)

        self.mycanvas = Canvas( width=300,height=300, bg='white')
        self.mycanvas.place(x=200,y=60)
      
        label = tk.Label(self.root, text='DIDO', font='Bold 20', bg='white')
        label.place(x=310,y=200)
        label = tk.Label(self.root, text='developed by Danai Bili, BSc.', font='Italics 8',bg='white')
        label.place(x=270,y=300)
        

    def getroot(self):
        return self.root

    def startup(self, fun):
        self.mycanvas = Canvas( width=300,height=300, bg='white')
        self.mycanvas.place(x=200,y=60)
        label = tk.Label(self.root, text='DIDO', font='Bold 15',bg='white')
        label.place(x=320,y=80)
        self.button = tk.Button(self.root,
                                    text="RUN",
                                    command=fun)
            
        self.button.pack()
        self.button.place(x=325,y=200)
        self.message = 'in'
 
    def dest_dir_choice(self):
       
        self.mycanvas = Canvas( width=300,height=300,bg='white')
        self.mycanvas.place(x=200,y=60)
        label = tk.Label(self.root, text='DIDO', font='Bold 15',bg='white')
        label.place(x=320,y=80)
        label = tk.Label(self.root, text='select destination folder', font='15',bg='white')
        label.place(x=270,y=180)
    

    def nirsdataload(self):
       
        self.mycanvas = Canvas( width=300,height=300,bg='white')
        self.mycanvas.place(x=200,y=60)
        label = tk.Label(self.root, text='DIDO', font='Bold 15',bg='white')
        label.place(x=320,y=80)
        label = tk.Label(self.root, text='select bNIRS file', font='15',bg='white')
        label.place(x=270,y=180)
    
    def dcsdataload(self):
        
        self.mycanvas = Canvas( width=300,height=300,bg='white')
        self.mycanvas.place(x=200,y=60)
        label = tk.Label(self.root, text='DIDO', font='Bold 15',bg='white')
        label.place(x=320,y=80)
        label = tk.Label(self.root, text='select DCS folder', font='15',bg='white')
        label.place(x=270,y=180)
        
    def nirsproc(self,nnirs):
    
        self.label = tk.Label(self.root, text='DIDO', font='Bold 15',bg='white')
        self.label.place(x=320,y=80)
        self.label = tk.Label(self.root, text="running nirs data processing", font='15',bg='white')
        self.label.place(x=240,y=180)
        self.root.update_idletasks()



    def dcsproc(self, ndcs):
        self.label = tk.Label(self.root, text='DIDO', font='Bold 15',bg='white')
        self.label.place(x=320,y=80)
        self.label = tk.Label(self.root, text="running dcs data processing", font='15',bg='white')
        self.label.place(x=240,y=180)
        #threading.Thread(target=self.bar(int( ndcs / 100))).start()
        self.root.update_idletasks()

    def accept_reject_buttons(self):

        self.button = tk.Button(self.root,
                                    text="accept",
                                    command=self.accept_option)
        self.button.pack()
        self.button.place(x=275,y=300)

        self.button = tk.Button(self.root,
                                    text="reject",
                                    command=self.reject_option)
        self.button.pack()
        self.button.place(x=375,y=300)

    def nirsdataaccept(self, accfun, rejfun):

        self.accept_option = accfun
        self.reject_option = rejfun
        
        self.mycanvas = Canvas( width=300,height=300,bg='white')
        self.mycanvas.place(x=200,y=60)
        label = tk.Label(self.root, text='DiDo', font='Bold 15',bg='white')
        label.place(x=320,y=80)
        label = tk.Label(self.root, text='Accept Dataset', font='15',bg='white')
        label.place(x=270,y=180)
        
        """img  = Image.open("code_development\\mus.png") 
        img = img.resize((160,100), Image.ANTIALIAS)
        photo=ImageTk.PhotoImage(img)
        lab=tk.Label(image=photo)  
        lab.place(x=200, y=200)"""
        #lab.lift()
        self.accept_reject_buttons()
            
    

        

    def dcsdataaccept(self, accfun, rejfun):
        self.accept_option = accfun
        self.reject_option = rejfun
        self.mycanvas = Canvas( width=300,height=300,bg='white')
        self.mycanvas.place(x=200,y=60)
        label = tk.Label(self.root, text='Dido', font='Bold 15',bg='white')
        label.place(x=320,y=80)
        label = tk.Label(self.root, text='Accept Dataset', font='15',bg='white')
        label.place(x=270,y=180)
        
        self.accept_reject_buttons()



    def mlproc(self):
        #self.message = p.runml()
        self.mycanvas = Canvas( width=300,height=300,bg='white')
        self.mycanvas.place(x=200,y=60)
        self.label = tk.Label(self.root, text='DIDO', font='Bold 15',bg='white')
        self.label.place(x=320,y=80)
        self.label = tk.Label(self.root, text="AI in action", font='15',bg='white')
        self.label.place(x=300,y=180)
        #self.bar(0.01)
        self.root.update_idletasks()

    def showres(self):
        #classif, conf = p.broadcastresult()
        classif = 'draft' ; conf = 'draft'; color='blue'
        if classif == 'mild':
            color = 'blue'
        if classif == 'moderate':
            color = 'orange'
        if classif == 'severe':
            color = 'red'
        self.mycanvas = Canvas( width=300,height=300,bg='white')
        self.mycanvas.place(x=200,y=60)
        label = tk.Label(self.root, text='DIDO', font='Bold 15',bg='white')
        label.place(x=320,y=80)
        label = tk.Label(self.root, text=f'diagnosis: {classif}', font='Bold 20', fg=color,bg='white')
        label.place(x=250,y=180)
        label = tk.Label(self.root, text=f'confidence: {conf}', font='Bold 16', bg='white')
        label.place(x=260,y=250)


    def bar(self,delay):
        def step():
            for i in range(int(100)):
                self.root.update_idletasks()
                pb1['value'] += 1 
                time.sleep(delay)
        pb1 = Progressbar(self.root, orient=HORIZONTAL, length=100, mode='indeterminate')
        pb1.pack(expand=True)
        pb1.place(x=295,y=250)
       #
        step()
        

