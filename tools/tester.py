def test_binary_sequences(file_template, test_class, output_file, *test_args):
    """
    Тестирует двоичные последовательности и сохраняет результаты в файл.

    :param file_template: Шаблон пути к файлам последовательностей (например, "../binary_sequences/sequence_{i}.txt").
    :param test_class: Класс теста (например, MDRNTest, MDRPTest и т.д.).
    :param output_file: Путь к выходному файлу для сохранения результатов теста.
    :param test_args: Аргументы, которые необходимы для инициализации теста (например, L, m, K и т.д.).
    """
    p_values = []
    alpha = 0.05

    # Загрузка последовательностей и проведение теста
    for i in range(100):
        filename = f"{file_template}_{i}.txt"
        with open(filename, "r") as sequence_file:
            sequence = sequence_file.read().strip()

            # Инициализация теста с передачей всех необходимых аргументов
            test_instance = test_class(sequence, *test_args)
            p_value = test_instance.run_test()

            p_values.append(p_value)
            print(f"Последовательность {i + 1}: p-значение = {p_value}")

    # Подсчет результатов
    passed_tests = sum(p > alpha for p in p_values)
    failed_tests = len(p_values) - passed_tests

    print(f"Passed: {passed_tests}")
    print(f"Not Passed: {failed_tests}")

    # Запись результатов в файл
    with open(output_file, "w") as result_file:
        result_file.write("p-values:\n")
        for i, p_value in enumerate(p_values, 1):
            if p_value <= alpha:
                result_file.write(f"Sequence {i}: p-value = {p_value} NOT PASSED\n")
            else:
                result_file.write(f"Sequence {i}: p-value = {p_value}\n")
        result_file.write(f"\nPassed: {passed_tests}\n")
        result_file.write(f"Not passed: {failed_tests}\n")

    print("Готово! Результаты теста и p-значения сохранены.")