import matplotlib.pyplot as plt

def get_pos_from_file(s1, s2, s3, s4):
    

    z = float(s1.split()[1])
    r = float(s2.split()[1])
    s4 = s4[1:-1]
    print(s4)
    re = float(s4.split(',')[0])
    im = float(s4.split(',')[1])

    return z, r, re, im, s3


H_rho_r_0_list_Im = list()
H_phi_r_0_list_Im = list()
H_z_r_0_list_Im = list()
E_rho_r_0_list_Im = list()
E_phi_r_0_list_Im = list()
E_z_r_0_list_Im = list()

H_rho_r_0_list_Re = list()
H_phi_r_0_list_Re = list()
H_z_r_0_list_Re = list()
E_rho_r_0_list_Re = list()
E_phi_r_0_list_Re = list()
E_z_r_0_list_Re = list()


H_rho_z_0_list_Im = list()
H_phi_z_0_list_Im = list()
H_z_z_0_list_Im = list()
E_rho_z_0_list_Im = list()
E_phi_z_0_list_Im = list()
E_z_z_0_list_Im = list()

H_rho_z_0_list_Re = list()
H_phi_z_0_list_Re = list()
H_z_z_0_list_Re = list()
E_rho_z_0_list_Re = list()
E_phi_z_0_list_Re = list()
E_z_z_0_list_Re = list()


for i in range(40):
    s1 = input()
    s2 = input()
    s3 = input()
    s4 = input()
    z, r, re, im, comp_name = get_pos_from_file(s1, s2, s3, s4)
    
    match comp_name:
        case "H_rho":
            if (r == 0.) :
                H_rho_r_0_list_Im.append((z, im))
                H_rho_r_0_list_Re.append((z, re))
            else:
                H_rho_z_0_list_Im.append((r, im))
                H_rho_z_0_list_Re.append((r, re))
        case "H_phi":
            if (r == 0.) :
                H_phi_r_0_list_Im.append((z, im))
                H_phi_r_0_list_Re.append((z, re))
            else:
                H_phi_z_0_list_Im.append((r, im))
                H_phi_z_0_list_Re.append((r, re))
        case "H_z":
            if (r == 0.):
                H_z_r_0_list_Im.append((z, im))
                H_z_r_0_list_Re.append((z, re))
            else:
                H_z_z_0_list_Im.append((r, im))
                H_z_z_0_list_Re.append((r, re))
        case "E_rho":
            if (r == 0.):
                E_rho_r_0_list_Im.append((z, im))
                E_rho_r_0_list_Re.append((z, re))
            else:
                E_rho_z_0_list_Im.append((r, im))
                E_rho_z_0_list_Re.append((r, re))
        case "E_z":
            if (r == 0.):
                E_z_r_0_list_Im.append((z, im))
                E_z_r_0_list_Re.append((z, re))
            else:
                E_z_z_0_list_Im.append((r, im))
                E_z_z_0_list_Re.append((r, re))
        case "E_phi":
            if (r == 0.):
                E_phi_r_0_list_Im.append((z, im))
                E_phi_r_0_list_Re.append((z, re))
            else:
                E_phi_z_0_list_Im.append((z, im))
                E_phi_z_0_list_Im.append((z, re))

dict_of_lists = dict()
dict_of_lists["H_rho_r_0_list_Im"] = H_rho_r_0_list_Im
dict_of_lists["H_phi_r_0_list_Im"] = H_phi_r_0_list_Im
dict_of_lists["H_z_r_0_list_Im"] = H_z_r_0_list_Im
dict_of_lists["E_rho_r_0_list_Im"] = E_rho_r_0_list_Im
dict_of_lists["E_phi_r_0_list_Im"] = E_phi_r_0_list_Im
dict_of_lists["E_z_r_0_list_Im"] = E_z_r_0_list_Im

dict_of_lists["H_rho_r_0_list_Re"] = H_rho_r_0_list_Re
dict_of_lists["H_phi_r_0_list_Re"] = H_phi_r_0_list_Re
dict_of_lists["H_z_r_0_list_Re"] = H_z_r_0_list_Re
dict_of_lists["E_rho_r_0_list_Re"] = E_rho_r_0_list_Re
dict_of_lists["E_phi_r_0_list_Re"] = E_phi_r_0_list_Re
dict_of_lists["E_z_r_0_list_Re"] = E_z_r_0_list_Re

dict_of_lists["H_rho_z_0_list_Im"] = H_rho_z_0_list_Im
dict_of_lists["H_phi_z_0_list_Im"]  = H_phi_z_0_list_Im
dict_of_lists["H_z_z_0_list_Im"] = H_z_z_0_list_Im
dict_of_lists["E_rho_z_0_list_Im"] = E_rho_z_0_list_Im
dict_of_lists["E_phi_z_0_list_Im"] = E_phi_z_0_list_Im
dict_of_lists["E_z_z_0_list_Im"] = E_z_z_0_list_Im

dict_of_lists["H_rho_z_0_list_Re"] = H_rho_z_0_list_Re
dict_of_lists["H_phi_z_0_list_Re"] = H_phi_z_0_list_Re
dict_of_lists["H_z_z_0_list_Re"] = H_z_z_0_list_Re
dict_of_lists["E_rho_z_0_list_Re"] = E_rho_z_0_list_Re
dict_of_lists["E_phi_z_0_list_Re"] = E_phi_z_0_list_Re
dict_of_lists["E_z_z_0_list_Re"] = E_z_z_0_list_Re

dict_of_names = dict()
dict_of_names["H_rho_r_0_list_Im"] = ("Im(H_rho) при r = 0", "z", "Im(H_rho)")
dict_of_names["H_phi_r_0_list_Im"] = ("Im(H_phi) при r = 0", "z", "Im(H_phi)")
dict_of_names["H_z_r_0_list_Im"] = ("Im(H_z) при r = 0", "z", "Im(H_z)")
dict_of_names["E_rho_r_0_list_Im"] = ("Im(E_rho) при r = 0", "z", "Im(E_rho)")
dict_of_names["E_phi_r_0_list_Im"] = ("Im(E_phi) при r = 0", "z", "Im(E_phi)")
dict_of_names["E_z_r_0_list_Im"] = ("Im(E_z) при r = 0", "z", "Im(E_z)")

dict_of_names["H_rho_r_0_list_Re"] = ("Re(H_rho) при r = 0", "z", "Re(H_rho)") 
dict_of_names["H_phi_r_0_list_Re"] = ("Re(H_phi) при r = 0", "z", "Re(H_phi)") 
dict_of_names["H_z_r_0_list_Re"] = ("Re(H_z) при r = 0", "z", "Re(H_z)") 
dict_of_names["E_rho_r_0_list_Re"] = ("Re(E_rho) при r = 0", "z", "Re(E_rho)") 
dict_of_names["E_phi_r_0_list_Re"] = ("Re(E_phi) при r = 0", "z", "Re(E_phi)") 
dict_of_names["E_z_r_0_list_Re"] = ("Re(E_z) при r = 0", "z", "Re(E_z)") 

dict_of_names["H_rho_z_0_list_Im"] = ("Im(H_rho) при z = z_0", "r", "Im(H_rho)") 
dict_of_names["H_phi_z_0_list_Im"]  = ("Im(H_phi) при z = z_0", "r", "Im(H_phi)")
dict_of_names["H_z_z_0_list_Im"] = ("Im(H_z) при z = z_0", "r", "Im(H_z)")
dict_of_names["E_rho_z_0_list_Im"] = ("Im(E_rho) при z = z_0", "r", "Im(E_rho)")
dict_of_names["E_phi_z_0_list_Im"] = ("Im(E_phi) при z = z_0", "r", "Im(E_phi)")
dict_of_names["E_z_z_0_list_Im"] = ("Im(E_z) при z = z_0", "r", "Im(E_z)")

dict_of_names["H_rho_z_0_list_Re"] = ("Re(H_rho) при z = z_0", "r", "Re(H_rho)")
dict_of_names["H_phi_z_0_list_Re"] = ("Re(H_phi) при z = z_0", "r", "Re(H_phi)")
dict_of_names["H_z_z_0_list_Re"] = ("Re(H_z) при z = z_0", "r", "Re(H_z)")
dict_of_names["E_rho_z_0_list_Re"] = ("Re(E_rho) при z = z_0", "r", "Re(E_rho)")
dict_of_names["E_phi_z_0_list_Re"] = ("Re(E_phi) при z = z_0", "r", "Re(E_phi)")
dict_of_names["E_z_z_0_list_Re"] = ("Re(E_z) при z = z_0", "r", "Re(E_z)")



for key in dict_of_lists.keys():
    x = list()
    y = list()
    preord = list()
    for arg, val in dict_of_lists[key]:
        preord.append((arg, val))
    preord = sorted(preord)
    for arg, val in preord:
        x.append(arg)
        y.append(val)
    
    if len(x) == 0:
        continue


    plt.plot(x, y, color='blue', marker='o')  # Линия с метками
    plt.title(dict_of_names[key][0])  # Заголовок графика
    plt.xlabel(dict_of_names[key][1])  # Метка оси X
    plt.ylabel(dict_of_names[key][2])  # Метка оси Y
    #plt.legend()  # Легенда
    plt.grid(True)  # Сетка
    plt.xscale("log")
    plt.show()  # Отображение графика
