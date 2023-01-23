
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
		enum class dataTypes : int { vector, matrix }; // { 0, 1 }
		std::vector<std::string> dataTypeNames = { "vector", "matrix"};
		std::vector<std::string> fileExtensions = { ".csv" };
		constexpr auto makeIndexingSet = [](int n) -> std::list<int> {
			std::list<int> ell(n);
			std::iota(ell.begin(), ell.end(), 0);
			return ell;
		};

		template <typename T> // cppreference.com/w/cpp/utility/to_underlying
		constexpr std::underlying_type_t<T> to_underlying(T v) noexcept
		{
			return static_cast<std::underlying_type_t<T>>(v);
		}

		constexpr auto isVector = [](int type) {
			bool binary_truth = type == to_underlying(dataTypes::vector);
			return binary_truth;
		};
		constexpr auto isMatrix = [](int type) {
			bool binary_truth = type == to_underlying(dataTypes::matrix);
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

		bool canRead(const std::string extension)
		{
			auto exts = makeIndexingSet(fileExtensions.size());
			for (auto index : exts)
				if (extension == fileExtensions[index])
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


		std::string getTerminalInput(const std::string prompt)
		{
			std::string input;
			std::cout << prompt << std::endl;
			std::getline(std::cin, input);
			return input;
		}

		std::optional<std::variant<vector, matrix>>
		getCSVDataFrom(const std::string feyell)
		{
			csv::CSVFormat format;
			format.header_row(0)
						.trim({' ', '\t'});
						//.variable_columns(true);
			csv::CSVReader reader(feyell);//, format);
			auto col_names = reader.get_col_names();
			auto indcs = makeIndexingSet(dataTypeNames.size());
			int type_of_data = -1; // I despise this method!

			for (auto col_name : col_names)
			{
				if (type_of_data >= 0)
					break;
				for (const auto indx : enum_range(dataTypes::vector, dataTypes::matrix))
				{
					int i = static_cast<int>(indx);
					if (col_name == dataTypeNames[i])
					{
						type_of_data = i;
						break;
					}
				}
			}
			if (type_of_data == -1)
				return std::nullopt;

			std::vector<std::vector<double>> rowsData;
			for (csv::CSVRow& row : reader)
			{
				std::vector<double> rowData;
    		for (csv::CSVField& field : row)
					rowData.push_back(field.get<double>());
				rowsData.push_back(rowData);
			}
			if (rowsData.empty())
				return std::nullopt;

			int num_of_rows = rowsData.size();
			int num_of_columns = rowsData.front().size();
			if (num_of_rows == 1)
			{
				std::vector<double> rowVector = rowsData.front();
				if (isVector(type_of_data))
				{
					vector v(rowVector);
					return v;
				} else if (isMatrix(type_of_data))
				{
					matrix m(rowVector);
					return m;
				}
			}
			else if (num_of_rows > 1 && isMatrix(type_of_data))
			{
				matrix m(rowsData);
				return m;
			}

			return std::nullopt;
		}

		// directory_iterator help from s.o. question:
		// 612097/how-can-i-get-the-list-of-files-in-a-directory-using-c-or-c
		std::optional<std::pair<std::vector<vector>, std::vector<matrix>>>
		getDataFrom(const std::string& directory)
		{
			std::vector<std::string> filenames;
    	for (const auto& entry : std::filesystem::directory_iterator(directory))
      {
				auto pathe = entry.path();
				if (!pathe.empty() && canRead(pathe.extension()))
					filenames.push_back(pathe);
			}
			if (filenames.empty())
				return std::nullopt;

			std::vector<vector> vees{};
			std::vector<matrix> emms{};

			for (auto feyell : filenames)
			{
				auto csvFileData = getCSVDataFrom(feyell);
				if (!csvFileData.has_value())
					break;
				auto v_or_m = std::move(csvFileData.value());
				const int i = v_or_m.index();
				if (isVector(i))
				{
					auto v = std::move(std::get<vector>(v_or_m));
					vees.push_back(std::move(v));
				} else if (isMatrix(i))
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