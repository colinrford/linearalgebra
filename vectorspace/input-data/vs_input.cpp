
#include "../matrix.hpp"
#include "../../submodules/csv-parser/single_include/csv.hpp"
#include <filesystem>
#include <list>
#include <numeric>
#include <ranges>
#include <sstream>
#include <type_traits>
#include <utility>
#include <variant>

namespace linalg
{
	namespace data_input
	{
		enum class data_types : int { vector, matrix }; // { 0, 1 }
		std::vector<std::string> data_type_names = { "vector", "matrix" };
		std::vector<std::string> file_extensions = { ".csv" };
		constexpr auto makeIndexingSet = [](int n) -> std::list<int> {
			std::list<int> ell(n);
			std::iota(ell.begin(), ell.end(), 0);
			return ell;
		};

		template <typename T> // en.cppreference.com/w/cpp/utility/to_underlying
		constexpr std::underlying_type_t<T> to_underlying(T v) noexcept
		{
			return static_cast<std::underlying_type_t<T>>(v);
		}

		constexpr auto is_vector = [](int type) {
			bool binary_truth = type == to_underlying(data_types::vector);
			return binary_truth;
		};
		constexpr auto is_matrix = [](int type) {
			bool binary_truth = type == to_underlying(data_types::matrix);
			return binary_truth;
		};

		/*
		auto constructVectorOrMatrix = [](int t1, int t2) {
			if (t1)
				do this;
			else if (t2)
				do that;
		};
		*/

		bool can_read(const std::string extension)
		{
			auto exts = makeIndexingSet(file_extensions.size());
			for (auto index : exts)
				if (extension == file_extensions[index])
					return true;
			return false;
		}

		constexpr inline auto enum_range = [](auto front, auto back)
		{
  		return std::views::iota(to_underlying(front),
															to_underlying(back) + 1)
       | std::views::transform([](auto e) { return decltype(front)(e); });
		};
		// enum_range lambda from stackoverflow post
		// /questions/69762598/what-are-commonly-used-ways-to-iterate-over-an-enum-class-in-c

		constexpr inline auto data_types_range = enum_range(data_types::vector,
																												data_types::matrix);

		std::string get_terminal_input(const std::string prompt)
		{
			std::string input;
			std::cout << prompt << std::endl;
			std::getline(std::cin, input);
			return input;
		}

		std::optional<std::variant<vector, matrix>>
		get_CSV_data_from(const std::string feyell)
		{
			csv::CSVFormat format;
			format.header_row(0)
						.trim({' ', '\t'});
						//.variable_columns(true);
			csv::CSVReader reader(feyell);//, format);
			auto col_names = reader.get_col_names();
			auto indcs = makeIndexingSet(data_type_names.size());
			int type_of_data = -1; // I despise this method!

			for (auto col_name : col_names)
			{
				if (type_of_data >= 0)
					break;
				for (const auto indx : data_types_range)
				{
					int i = static_cast<int>(indx);
					if (col_name == data_type_names[i])
					{
						type_of_data = i;
						break;
					}
				}
			}
			if (type_of_data == -1)
				return std::nullopt;

			std::vector<std::vector<double>> rows_data;
			for (csv::CSVRow& row : reader)
			{
				std::vector<double> row_data;
    		for (csv::CSVField& field : row)
					row_data.push_back(field.get<double>());
				rows_data.push_back(row_data);
			}
			if (rows_data.empty())
				return std::nullopt;

			int num_of_rows = rows_data.size();
			int num_of_columns = rows_data.front().size();
			if (num_of_rows == 1)
			{
				std::vector<double> row_vector = rows_data.front();
				if (is_vector(type_of_data))
				{
					vector v(row_vector);
					return v;
				} else if (is_matrix(type_of_data))
				{
					matrix m(row_vector);
					return m;
				}
			}
			else if (num_of_rows > 1 && is_matrix(type_of_data))
			{
				matrix m(rows_data);
				return m;
			}

			return std::nullopt;
		}

		// directory_iterator help from s.o. question:
		// 612097/how-can-i-get-the-list-of-files-in-a-directory-using-c-or-c
		std::optional<std::pair<std::vector<vector>, std::vector<matrix>>>
		get_data_from(const std::string& directory)
		{
			std::vector<std::string> filenames;
    	for (const auto& entry : std::filesystem::directory_iterator(directory))
      {
				auto pathe = entry.path();
				if (!pathe.empty() && can_read(pathe.extension()))
					filenames.push_back(pathe);
			}
			if (filenames.empty())
				return std::nullopt;

			std::vector<vector> vees{};
			std::vector<matrix> emms{};

			for (auto feyell : filenames)
			{
				auto csvFileData = get_CSV_data_from(feyell);
				if (!csvFileData.has_value())
					break;
				auto v_or_m = std::move(csvFileData.value());
				const int i = v_or_m.index();
				if (is_vector(i))
				{
					auto v = std::move(std::get<vector>(v_or_m));
					vees.push_back(std::move(v));
				} else if (is_matrix(i))
				{
					auto m = std::move(std::get<matrix>(v_or_m));
					emms.push_back(std::move(m));
				}
			}

			if (vees.empty() && emms.empty())
				return std::nullopt;

			auto pear = std::make_pair(std::move(vees), std::move(emms));
			return pear;
		}

	}
}
