def lfsr_to_file(seed, taps, num_bits, output_file):
    """
    Генерация последовательности LFSR и запись в файл.

    :param seed: Начальное состояние (список битов, например, [1, 0, 0, 1]).
    :param taps: Позиции обратной связи (список индексов, например, [0, 2]).
    :param num_bits: Количество битов в выходной последовательности.
    :param output_file: Имя файла для записи выходной последовательности.
    """
    # Проверяем, чтобы длина seed была достаточной
    state = seed[:]  # Копируем начальное состояние
    n = len(state)  # Длина регистра

    # Проверяем, что все индексы из taps находятся в пределах регистра
    if any(tap >= n or tap < 0 for tap in taps):
        raise ValueError(f"Некорректные индексы в taps: {taps}. Длина регистра: {n}.")

    with open(output_file, "wb") as f:
        byte = 0  # Текущий байт
        bit_count = 0  # Счётчик битов в байте

        for _ in range(num_bits):
            # Берём самый правый бит как выходной
            output_bit = state[-1]

            # Добавляем бит в текущий байт
            byte = (byte << 1) | output_bit
            bit_count += 1

            # Если накопили 8 бит, записываем байт в файл
            if bit_count == 8:
                f.write(byte.to_bytes(1, byteorder="big"))
                byte = 0
                bit_count = 0

            # Рассчитываем новый бит с помощью XOR
            new_bit = 0
            for tap in taps:
                new_bit ^= state[-(tap + 1)]  # XOR всех указанных битов обратной связи

            # Сдвигаем регистр влево и добавляем новый бит
            state = [new_bit] + state[:-1]

        # Если остались незаписанные биты, записываем их (дополняем нулями до байта)
        if bit_count > 0:
            byte = byte << (8 - bit_count)  # Дополняем оставшиеся биты нулями
            f.write(byte.to_bytes(1, byteorder="big"))


# Пример использования
if __name__ == "__main__":
    # Параметры LFSR
    degree = 30  # Степень LFSR
    seed = [1] * degree
    taps = [0, 1, 9, 29]  # Обратная связь для полинома x^30 + x^29 + x^9 + x^0
    num_bits = 100 * 1024 * 1024 * 8
    output_file = "../sequences/lfsr(x^30+x^29+x^9+1).bin"  # Имя выходного файла

    # Генерация последовательности
    print(f"Генерация последовательности длиной {num_bits} бит...")
    lfsr_to_file(seed, taps, num_bits, output_file)
    print(f"Последовательность записана в файл: {output_file}")