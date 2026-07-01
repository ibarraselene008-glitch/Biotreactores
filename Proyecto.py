import tkinter as tk
from tkinter import messagebox, ttk
import matplotlib.pyplot as plt

DICCIONARIO_MICROORGANISMOS = {
    "Geobacillus stearothermophilus (Esporas)": {"Ea": 282343, "A": 4.93e37},
    "Bacillus subtilis (Esporas)": {"Ea": 243509, "A": 9.5e32},
    "Gluconobacter oxydans (Célula vegetativa)": {"Ea": 182004, "A": 1.8e24},
    "Escherichia coli (Célula vegetativa)": {"Ea": 172381, "A": 1.2e25},
    "Otro (Ingresar manualmente)": {"Ea": 0.0, "A": 0.0},
}

def actualizar_constantes_biol(event):
    """Actualiza las casillas de texto de Ea y A según el microorganismo seleccionado."""
    seleccionado = combo_microorganismo.get()
    datos = DICCIONARIO_MICROORGANISMOS[seleccionado]

    entry_Ea.config(state="normal")
    entry_Aarr.config(state="normal")

    entry_Ea.delete(0, tk.END)
    entry_Ea.insert(0, str(datos["Ea"]))

    entry_Aarr.delete(0, tk.END)
    entry_Aarr.insert(0, f"{datos['A']:.2e}")

    if seleccionado != "Otro (Ingresar manualmente)":
        entry_Ea.config(state="disabled")
        entry_Aarr.config(state="disabled")

def simular_esterilizacion():
    try:
        n0 = float(entry_n0.get())
        Vop = float(entry_Vop.get())
        N_final_obj = float(entry_Nfinal.get())

        N0 = n0 * Vop

        if N0 <= 0 or N_final_obj <= 0:
            raise ValueError

        nabla_requerido = math.log(N0 / N_final_obj)

        Ea = float(entry_Ea.get())
        A_arr = float(entry_Aarr.get())
        R = 8.314

        if Ea <= 0 or A_arr <= 0:
            raise ValueError

        Tmdeseada = float(entry_Tmdeseada.get()) + 273.15
        Tenfriamiento = float(entry_Tenfriamiento.get()) + 273.15
        T0 = float(entry_T0.get()) + 273.15
        Th = float(entry_Th.get()) + 273.15
        Tco = float(entry_Tco.get()) + 273.15

        # VALIDACIONES FÍSICAS
        if T0 >= Tmdeseada:
            messagebox.showerror(
                "Error de Temperaturas",
                "La temperatura inicial debe ser menor que la temperatura de mantenimiento."
            )
            return

        if Tmdeseada >= Th:
            messagebox.showerror(
                "Error de Temperaturas",
                "La temperatura de mantenimiento debe ser menor que la temperatura del vapor."
            )
            return

        if Tenfriamiento <= Tco:
            messagebox.showerror(
                "Error de Temperaturas",
                "La temperatura final de enfriamiento debe ser mayor que la temperatura del agua de enfriamiento."
                            )
            return

        U_cal = float(entry_Ucal.get())
        U_enf = float(entry_Uenf.get())
        Area = float(entry_Area.get())
        Masa = float(entry_Masa.get())
        Cp = float(entry_Cp.get())
        W = float(entry_W.get())
        Cp_agua = float(entry_Cpagua.get())

        if U_cal <= 0 or U_enf <= 0 or Area <= 0 or Masa <= 0 or Cp <= 0 or W <= 0 or Cp_agua <= 0:
            messagebox.showerror(
                "Error de Parámetros",
                "Todos los parámetros deben ser positivos."
            )
            return

        alpha = ((U_cal * Area) / (Masa * Cp))*60  # min^-1
        beta = (T0 - Th) / Th
        upsilon = ((W * Cp_agua) / (Masa * Cp) * (
            1 - math.exp(-(U_enf * Area) / (W * Cp_agua))
        ))*60  # min^-1
        T_inicio_enf = Tmdeseada
        lambda_param = (T_inicio_enf - Tco) / Tco

        dt = 0.001  # min

        tiempos_grafica = []
        temperaturas_grafica = []
        nablas_grafica = []

        t = 0.0
        T_actual = T0
        nabla_C = 0.0

        tiempos_grafica.append(t)
        temperaturas_grafica.append(T_actual - 273.15)
        nablas_grafica.append(nabla_C)

        while T_actual < Tmdeseada :
            t += dt
            T_actual = Th * (1.0 + beta * math.exp(-alpha * t))

    
            if T_actual > (Tmdeseada):
                T_actual = Tmdeseada

            k = A_arr * math.exp(-Ea / (R * T_actual))
            nabla_C += k * dt

            if abs(t % 1.0) < dt:
                tiempos_grafica.append(round(t))
                temperaturas_grafica.append(T_actual - 273.15)
                nablas_grafica.append(nabla_C)

        t_calentamiento = t
        if T_actual < Tmdeseada:
            messagebox.showwarning(
                "Alerta de Calentamiento",
                "No se alcanzó la temperatura de mantenimiento dentro del tiempo simulado."
            )
            return

        t_e = 0.0
        T_enf_actual = T_inicio_enf
        nabla_E = 0.0

        while T_enf_actual > Tenfriamiento:
            t_e += dt
            T_enf_actual = Tco * (1.0 + lambda_param * math.exp(-upsilon * t_e))
            k = A_arr * math.exp(-Ea / (R * T_enf_actual))
            nabla_E += k * dt

        nabla_M = (
            nabla_requerido - nabla_C - nabla_E
            if (nabla_requerido - nabla_C - nabla_E) > 0
            else 0.0
        )
        k_man = A_arr * math.exp(-Ea / (R * (Tmdeseada)))
        t_mantenimiento = nabla_M / k_man if k_man > 0 else 0.0

        t_global = t_calentamiento
        nabla_total = nabla_C

        minutos_mantenimiento = int(math.ceil(t_mantenimiento/dt))
        for _ in range(minutos_mantenimiento):
            t_global += dt
            nabla_total += k_man * dt
            if abs(t_global % 1.0) < dt:
                tiempos_grafica.append(round(t_global))
                temperaturas_grafica.append(Tmdeseada - 273.15)
                nablas_grafica.append(nabla_total)

        t_e_real = 0.0
        T_enf_actual = T_inicio_enf
        while T_enf_actual > (Tenfriamiento) and t_e_real < t_e:
            t_e_real += dt
            t_global += dt
            T_enf_actual = Tco * (1.0 + lambda_param * math.exp(-upsilon * t_e_real))
            k = A_arr * math.exp(-Ea / (R * T_enf_actual))
            nabla_total += k * dt

            if abs(t_e_real % 1.0) < dt:
                tiempos_grafica.append(round(t_global))
                temperaturas_grafica.append(T_enf_actual - 273.15)
                nablas_grafica.append(nabla_total)

        N_final_real = N0 * math.exp(-nabla_total)

  
        porcentaje_C = (nabla_C / nabla_total) * 100 if nabla_total > 0 else 0
        porcentaje_M = (nabla_M / nabla_total) * 100 if nabla_total > 0 else 0
        porcentaje_E = (nabla_E / nabla_total) * 100 if nabla_total > 0 else 0

        mostrar_tabla(tiempos_grafica, temperaturas_grafica, nablas_grafica)

        messagebox.showinfo(
            "Reporte de Diseño Cinético",
            f"REPORTE ACADÉMICO DE ESTERILIZACIÓN:\n\n"
            f"• Microorganismo evaluado: {combo_microorganismo.get()}\n"
            f"  [Ea = {Ea:.0f} J/mol | A = {A_arr:.2e} min⁻¹]\n"
            f"-----------------------------------------\n"
            f"• Carga Inicial Total (N₀): {N0:.2e} UFC\n"
            f"• Población Superviviente Real (N): {N_final_real:.2e} UFC\n"
            f"-----------------------------------------\n"
            f"• Nabla Requerido (∇ de diseño): {nabla_requerido:.4f}\n"
            f"• Nabla Calculado (∇ total neto): {nabla_total:.4f}\n\n"
            f"DISTRIBUCIÓN DE LETALIDAD CINÉTICA:\n"
            f"• ∇ Calentamiento (∇C): {nabla_C:.4f} ({porcentaje_C:.1f}% del total)\n"
            f"• ∇ Mantenimiento (∇M): {nabla_M:.4f} ({porcentaje_M:.1f}% del total)\n"
            f"• ∇ Enfriamiento (∇E): {nabla_E:.4f} ({porcentaje_E:.1f}% del total)\n"
            f"-----------------------------------------\n"
            f"• Tiempo Calentamiento (tc): {t_calentamiento:.2f} min\n"
            f"• Tiempo Mantenimiento (tm): {t_mantenimiento:.2f} min\n"
            f"• Tiempo Enfriamiento (te): {t_e:.2f} min\n"
            f"• Tiempo Total de Operación: {t_global:.2f} min",
        )

        # Gráfica del perfil
        plt.figure(figsize=(9, 5))
        plt.plot(
            tiempos_grafica,
            temperaturas_grafica,
            color="#008080",
            linewidth=2.5,
            label="Temperatura del Caldo",
        )
        plt.axhline(
            y=Tmdeseada - 273.15, color="r", linestyle="--", label="Temperatura de Diseño (°C)"
        )
        plt.title(
            f"Perfil Analítico de Esterilización Térmica - Criterio ∇\n[Microorganismo: {combo_microorganismo.get()}]",
            fontsize=11,
            fontweight="bold",
        )
        plt.xlabel("Tiempo de Operación (min)", fontsize=10)
        plt.ylabel("Temperatura Interior del Biorreactor (°C)", fontsize=10)
        plt.grid(True, linestyle=":")
        plt.legend(loc="lower right")
        plt.show()

    except ValueError:
        messagebox.showerror(
            "Error de Parámetros",
            "Por favor, introduzca valores numéricos válidos mayores a cero.",
        )

def mostrar_tabla(tiempos, temperaturas, nablas):
    ventana_tabla = tk.Toplevel()
    ventana_tabla.title("Integración del Criterio Nabla")
    ventana_tabla.geometry("460x400")

    frame_tabla = ttk.Frame(ventana_tabla, padding="10")
    frame_tabla.pack(fill=tk.BOTH, expand=True)

    columnas = ("tiempo", "temp", "nabla")
    tabla = ttk.Treeview(frame_tabla, columns=columnas, show="headings")
    tabla.heading("tiempo", text="Tiempo aprox (min)")
    tabla.heading("temp", text="Temperatura (°C)")
    tabla.heading("nabla", text="∇ Acumulado Calculado")

    tabla.column("tiempo", anchor=tk.CENTER, width=110)
    tabla.column("temp", anchor=tk.CENTER, width=120)
    tabla.column("nabla", anchor=tk.CENTER, width=160)

    scrollbar = ttk.Scrollbar(frame_tabla, orient=tk.VERTICAL, command=tabla.yview)
    tabla.configure(yscrollcommand=scrollbar.set)

    for t, temp, nab in zip(tiempos, temperaturas, nablas):
        tabla.insert("", tk.END, values=(t, f"{temp:.2f}", f"{nab:.4f}"))

    tabla.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
    scrollbar.pack(side=tk.RIGHT, fill=tk.Y)


# INTERFAZ GRÁFICA 

root = tk.Tk()
root.title("Simulador de Esterilización Térmica por Lotes - UPIBI IPN 2.0")
root.geometry("560x640")

main_frame = ttk.Frame(root, padding="15")
main_frame.pack(fill=tk.BOTH, expand=True)

# SECCIÓN 1: MICROBIOLÓGICA
ttk.Label(
    main_frame, text="1.PARÁMETROS CINÉTICOS", font=("Arial", 10, "bold")
).grid(row=0, column=0, columnspan=2, sticky=tk.W, pady=(5, 5))

ttk.Label(main_frame, text="Microorganismo de Referencia:").grid(
    row=1, column=0, sticky=tk.W, pady=2
)
combo_microorganismo = ttk.Combobox(
    main_frame, values=list(DICCIONARIO_MICROORGANISMOS.keys()), state="readonly", width=35
)
combo_microorganismo.grid(row=1, column=1, sticky=tk.E, pady=2)
combo_microorganismo.bind("<<ComboboxSelected>>", actualizar_constantes_biol)

ttk.Label(main_frame, text="Energía de Activación (Ea, J/mol):").grid(
    row=2, column=0, sticky=tk.W, pady=2
)
entry_Ea = ttk.Entry(main_frame, width=15)
entry_Ea.grid(row=2, column=1, sticky=tk.E, pady=2)

ttk.Label(main_frame, text="Factor Preexponencial A (min^-1):").grid(
    row=3, column=0, sticky=tk.W, pady=2
)
entry_Aarr = ttk.Entry(main_frame, width=15)
entry_Aarr.grid(row=3, column=1, sticky=tk.E, pady=2)

ttk.Label(main_frame, text="Concentración Inicial (n0, UFC/mL):").grid(
    row=4, column=0, sticky=tk.W, pady=2
)
entry_n0 = ttk.Entry(main_frame, width=15)
entry_n0.insert(0, "1e6")
entry_n0.grid(row=4, column=1, sticky=tk.E, pady=2)

ttk.Label(main_frame, text="Volumen de Operación (Vop, mL):").grid(
    row=5, column=0, sticky=tk.W, pady=2
)
entry_Vop = ttk.Entry(main_frame, width=15)
entry_Vop.insert(0, "1000000")
entry_Vop.grid(row=5, column=1, sticky=tk.E, pady=2)

ttk.Label(main_frame, text="Población Superviviente Deseada (N, UFC):").grid(
    row=6, column=0, sticky=tk.W, pady=2
)
entry_Nfinal = ttk.Entry(main_frame, width=15)
entry_Nfinal.insert(0, "0.001")
entry_Nfinal.grid(row=6, column=1, sticky=tk.E, pady=2)

combo_microorganismo.current(0)
actualizar_constantes_biol(None)

ttk.Separator(main_frame, orient="horizontal").grid(
    row=7, column=0, columnspan=2, sticky="ew", pady=10
)

# SECCIÓN 2: INGENIERÍA
ttk.Label(
    main_frame, text="2. CAPACIDAD TÉRMICA Y DIMENSIONES DEL REACTOR", font=("Arial", 10, "bold")
).grid(row=8, column=0, columnspan=2, sticky=tk.W, pady=(0, 5))

campos_ingenieria = [
    ("Temp. deseada de mantenimiento (°C):", "121", "entry_Tmdeseada"),
    ("Temp. deseada de enfriamiento (°C):", "40", "entry_Tenfriamiento"),
    ("Temp. Inicial del Medio (°C):", "20", "entry_T0"),
    ("Temp. Vapor de Servicio (Th en °C):", "134", "entry_Th"),
    ("Temp. Agua Enfriamiento (Tco en °C):", "20", "entry_Tco"),
    ("Coef. U Calentamiento (W/m² °K):", "800", "entry_Ucal"),
    ("Coef. U Enfriamiento (W/m² °K):", "500", "entry_Uenf"),
    ("Área de Transferencia (m²):", "4.18", "entry_Area"),
    ("Masa del Medio (Kg):", "1134", "entry_Masa"),
    ("Cp del Medio (J/kg °K):", "4180", "entry_Cp"),
    ("Flujo de Agua Enfriamiento (kg/min):", "54.4", "entry_W"),
    ("Cp Agua Enfriamiento (J/kg °K):", "4180", "entry_Cpagua"),
]

for idx, (label_text, default_val, var_name) in enumerate(campos_ingenieria):
    current_row = idx + 9
    ttk.Label(main_frame, text=label_text).grid(row=current_row, column=0, sticky=tk.W, pady=2)
    entry = ttk.Entry(main_frame, width=15)
    entry.insert(0, default_val)
    entry.grid(row=current_row, column=1, sticky=tk.E, pady=2)
    globals()[var_name] = entry

btn_simular = ttk.Button(
    main_frame, text="Diseñar Proceso Térmico (Calcular ∇)", command=simular_esterilizacion
)
btn_simular.grid(
    row=len(campos_ingenieria) + 9,
    column=0,
    columnspan=2,
    pady=15,
    ipadx=15,
    ipady=6,
)

root.mainloop()


