#ifndef _MAGNETIC_FIELD_LUKE_H
#define _MAGNETIC_FIELD_LUKE_H

#include <string>
#include <softlib/MagneticField/MagneticFieldNumeric2D.h>

class MagneticFieldLUKE : public MagneticFieldNumeric2D {
	public:
		MagneticFieldLUKE(const std::string&);
		MagneticFieldLUKE(const std::string&, enum sfile_type);
		MagneticFieldLUKE(const std::string&, const std::string&);
		MagneticFieldLUKE(const std::string&, enum sfile_type, const std::string&);
		MagneticFieldLUKE(const std::string&, enum sfile_type, const std::string&, enum sfile_type);

		virtual void Load(const std::string&) override;
		virtual void Load(const std::string&, enum sfile_type) override;
		virtual void Load(const std::string&, const std::string&);
		virtual void Load(const std::string&, enum sfile_type, const std::string&);
		virtual void Load(const std::string&, enum sfile_type, const std::string&, enum sfile_type);

		virtual void LoadWall(const std::string&);
		virtual void LoadWall(const std::string&, enum sfile_type);
};

#endif/*_MAGNETIC_FIELD_LUKE_H*/
