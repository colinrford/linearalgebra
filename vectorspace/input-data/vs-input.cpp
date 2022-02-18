
#include "../matrix.h"
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
		enum class dataTypes { vector, matrix }; // { 0, 1 }
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

		/*
		auto constructVectorOrMatrix = [](int t1, int t2) {
			if (t1)
				do this;
			else if (t2)
				do that;
		};
		*/

		bool canRead(std::string extension)
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

		std::optional<std::variant<Vector, Matrix>>
		getCSVDataFrom(const std::string feyell)
		{
			csv::CSVFormat format;
			format.header_row(0)
						.trim({' ', '\t'});
						//.variable_columns(true);
			csv::CSVReader reader(feyell, format);
			auto col_names = reader.get_col_names();
			auto indcs = makeIndexingSet(dataTypeNames.size());
			int type_of_data = -1;

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
				std::cout << "is row?" << std::endl;
				std::vector<double> rowData;
    		for (csv::CSVField& field : row)
        {
					std::cout << "field.get<>() = " << field.get<>() << std::endl;
					rowData.push_back(field.get<double>());
				}
				std::cout << "numba tzu" << std::endl;
				rowsData.push_back(rowData);
			}
			if (rowsData.empty())
				return std::nullopt;

			int num_of_rows = rowsData.size();
			int num_of_columns = rowsData.front().size();
			if (num_of_rows == 1)
			{
				std::vector<double> rowVector = rowsData.front();
				if (type_of_data == 0)
				{
					Vector v(rowVector);
					return std::move(v);
				} else if (type_of_data == 1)
				{
					Matrix m(rowVector);
					return std::move(m);
				}
			}
			else if (num_of_rows > 1 && type_of_data == 1)
			{
				Matrix m(rowsData);
				return std::move(m);
			}

			return std::nullopt;
		}

		// directory_iterator help from s.o. question:
		// 612097/how-can-i-get-the-list-of-files-in-a-directory-using-c-or-c
		std::optional<std::pair<std::vector<Vector>, std::vector<Matrix>>>
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

			std::vector<Vector> vees{};
			std::vector<Matrix> emms{};

			for (auto feyell : filenames)
			{
				auto csvFileData = getCSVDataFrom(feyell);
				if (!csvFileData.has_value())
					break;
				auto v_or_m = std::move(csvFileData.value());
				const int i = v_or_m.index();
				if (i == static_cast<int>(dataTypes::vector))
				{
					auto v = std::move(std::get<Vector>(v_or_m));
					vees.push_back(std::move(v));
				}
				else if (i == static_cast<int>(dataTypes::matrix))
				{
					auto m = std::move(std::get<Matrix>(v_or_m));
					emms.push_back(std::move(m));
				}
			}

			if (vees.empty() && emms.empty())
				return std::nullopt;

			auto pear = std::make_pair(std::move(vees), std::move(emms));
			return std::move(pear);
		}


		void wip()
		{
			auto realm = getTerminalInput("directory of data files pls");
			auto possibilities = getDataFrom(realm);
			if (!possibilities.has_value())
				return;
			auto dataPear = std::move(possibilities.value());
			auto vectorData = std::move(std::get<std::vector<Vector>>(dataPear));
			auto matrixData = std::move(std::get<std::vector<Matrix>>(dataPear));

		}

	}
}
