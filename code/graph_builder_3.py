import re
import matplotlib.pyplot as plt
import numpy as np

def parse_data(text):
    lines = text.split('\n')
    pattern = re.compile(
        r'z:\s*(?P<z>[\d.-]+)\s*R:\s*(?P<R>[\d.-]+)\s*E_phi\s*\{\s*(?P<re>[\d.-]+),\s*(?P<im>[\d.-]+)\s*}'
    )
    
    # Извлекаем имя графика (3-я строка и повторяется каждые 3+4n строк)
    graph_names = [line.strip() for i, line in enumerate(lines) if i % 7 == 2 and line.strip()]
    graph_name = graph_names[0] if graph_names else "E_phi Components"
    
    z_values = []
    re_values = []
    im_values = []
    
    for match in pattern.finditer(text):
        z_values.append(float(match.group('z')))
        re_values.append(float(match.group('re')))
        im_values.append(float(match.group('im')))
    
    return z_values, re_values, im_values, graph_name

def plot_data(z_values, re_values, im_values, graph_name):
    plt.figure(figsize=(14, 6))
    plt.suptitle(graph_name, fontsize=14, y=1.02)
    
    # График действительной части
    plt.subplot(1, 2, 1)
    plt.semilogx(z_values, re_values, 'b-', label='Re(E_phi)')
    plt.semilogx(z_values, re_values, 'bo', markersize=5, alpha=0.7, label='Точки Re(E_phi)')
    plt.xlabel('z (log scale)')
    plt.ylabel('Re(E_phi)')
    plt.title('Действительная часть')
    plt.grid(True, which="both", ls="-")
    plt.legend()
    
    # График мнимой части
    plt.subplot(1, 2, 2)
    plt.semilogx(z_values, im_values, 'r-', label='Im(E_phi)')
    plt.semilogx(z_values, im_values, 'ro', markersize=5, alpha=0.7, label='Точки Im(E_phi)')
    plt.xlabel('z (log scale)')
    plt.ylabel('Im(E_phi)')
    plt.title('Мнимая часть')
    plt.grid(True, which="both", ls="-")
    plt.legend()
    
    plt.tight_layout()
    plt.show()

def main():
    # Чтение данных из файла или стандартного ввода
    print("Введите данные (Ctrl+D или Ctrl+Z для завершения ввода):")
    data = []
    while True:
        try:
            line = input()
            data.append(line)
        except EOFError:
            break
    
    text = '\n'.join(data)
    z_values, re_values, im_values, graph_name = parse_data(text)
    
    if not z_values:
        print("Не найдено данных для построения графиков.")
        return
    
    # Проверка и обработка отрицательных или нулевых значений z
    z_values = np.array(z_values)
    if np.any(z_values <= 0):
        print("Предупреждение: найдены неположительные значения z, заменяю на малые положительные.")
        min_positive_z = np.min(z_values[z_values > 0])
        z_values[z_values <= 0] = min_positive_z / 100
    
    plot_data(z_values, re_values, im_values, graph_name)

if __name__ == "__main__":
    main()