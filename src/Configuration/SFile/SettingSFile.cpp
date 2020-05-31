/**
 * Implementation of the SFile Setting.
 */

#include <limits>
#include <vector>
#include <softlib/Configuration/ConfigurationSFile.h>

using namespace std;


/**
 * Destructor.
 */
SettingSFile::~SettingSFile() {
    DeleteCurrent();
}

void SettingSFile::DeleteCurrent() {
    switch (this->dtype) {
        case DataType::BOOL:   delete this->data.b; break;
        case DataType::INT32:  delete this->data.i32; break;
        case DataType::INT64:  delete this->data.i64; break;
        case DataType::UINT32: delete this->data.u32; break;
        case DataType::UINT64: delete this->data.u64; break;
        case DataType::REAL:   delete this->data.r; break;
        case DataType::STRING: delete this->data.s; break;

        default: break;
    }

    this->dtype = DataType::UNSPECIFIED;
}

/**
 * Overwrite the values of this setting with the
 * values in the given setting.
 *
 * s: Setting, which's values are to overwrite the values
 *    of this setting.
 */
void SettingSFile::OverwriteValues(const Setting *s) {
    if (typeid(s) == typeid(SettingSFile*)) {
        const SettingSFile *ss = (SettingSFile*)s;

        switch (ss->dtype) {
            case DataType::BOOL:   this->SetValueBool(new vector<bool>(*(ss->data.b))); break;
            case DataType::INT32:  this->SetValueInt32(new vector<int32_t>(*(ss->data.i32))); break;
            case DataType::INT64:  this->SetValueInt64(new vector<int64_t>(*(ss->data.i64))); break;
            case DataType::UINT32: this->SetValueUInt32(new vector<uint32_t>(*(ss->data.u32))); break;
            case DataType::UINT64: this->SetValueUInt64(new vector<uint64_t>(*(ss->data.u64))); break;
            case DataType::REAL:   this->SetValueReal(new vector<slibreal_t>(*(ss->data.r))); break;
            case DataType::STRING: this->SetValueString(new vector<string>(*(ss->data.s))); break;

            default:
                throw ConfigurationException("The data in setting '%s' is of unrecognized type.", this->name.c_str());
        }
    } else {
        if (s->IsBool()) {
            vector<bool> *v = new vector<bool>(s->GetNumberOfValues());
            for (size_t i = 0; i < s->GetNumberOfValues(); i++)
                v->push_back(s->GetBool(i));

            this->SetValueBool(v);
        } else if (s->IsInteger32()) {
            vector<int32_t> *v = new vector<int32_t>(s->GetNumberOfValues());
            for (size_t i = 0; i < s->GetNumberOfValues(); i++)
                v->push_back(s->GetInteger32(i));

            this->SetValueInt32(v);
        } else if (s->IsInteger64()) {
            vector<int64_t> *v = new vector<int64_t>(s->GetNumberOfValues());
            for (size_t i = 0; i < s->GetNumberOfValues(); i++)
                v->push_back(s->GetInteger64(i));

            this->SetValueInt64(v);
        } else if (s->IsUnsignedInteger32()) {
            vector<uint32_t> *v = new vector<uint32_t>(s->GetNumberOfValues());
            for (size_t i = 0; i < s->GetNumberOfValues(); i++)
                v->push_back(s->GetUnsignedInteger32(i));

            this->SetValueUInt32(v);
        } else if (s->IsUnsignedInteger64()) {
            vector<uint64_t> *v = new vector<uint64_t>(s->GetNumberOfValues());
            for (size_t i = 0; i < s->GetNumberOfValues(); i++)
                v->push_back(s->GetUnsignedInteger64(i));

            this->SetValueUInt64(v);
        } else if (s->IsNumericVector()) {
            vector<slibreal_t> *v = new vector<slibreal_t>(s->GetNumberOfValues());
            for (size_t i = 0; i < s->GetNumberOfValues(); i++)
                v->push_back(s->GetScalar(i));

            this->SetValueReal(v);
        } else { // String
            vector<string> *v = new vector<string>(s->GetNumberOfValues());
            for (size_t i = 0; i < s->GetNumberOfValues(); i++)
                v->push_back(s->GetString());

            this->SetValueString(v);
        }
    }
}

/**
 * Make a copy of this setting.
 */
Setting *SettingSFile::Copy() {
    switch (this->dtype) {
        case DataType::BOOL:   return new SettingSFile(new vector<bool>(*(this->data.b)));
        case DataType::INT32:  return new SettingSFile(new vector<int32_t>(*(this->data.i32)));
        case DataType::INT64:  return new SettingSFile(new vector<int64_t>(*(this->data.i64)));
        case DataType::UINT32: return new SettingSFile(new vector<uint32_t>(*(this->data.u32)));
        case DataType::UINT64: return new SettingSFile(new vector<uint64_t>(*(this->data.u64)));

        default:
            throw ConfigurationException("The data in setting '%s' is of unrecognized type.", this->name.c_str());
    }
}

/************************
 * GETTERS              *
 ************************/
bool SettingSFile::GetBool(unsigned int index) const {
    if (IsBool())
        return this->data.b->at(index);
    else
        throw ConfigurationException("The data in setting '%s' is not of boolean type.", this->name.c_str());
}
int64_t SettingSFile::GetInteger(unsigned int index) const  {
    switch (this->dtype) {
        case DataType::INT32:  return this->data.i32->at(index);
        case DataType::INT64:  return this->data.i64->at(index);
        case DataType::UINT32: return this->data.u32->at(index);
        case DataType::UINT64: return this->data.u64->at(index);
        case DataType::REAL:   return this->data.r->at(index);

        default:
            throw ConfigurationException("The data in setting '%s' is not of integer type.", this->name.c_str());
    }
}
int32_t SettingSFile::GetInteger32(unsigned int index) const { return (int32_t)GetInteger(index); }
uint32_t SettingSFile::GetUnsignedInteger32(unsigned int index) const {return (uint32_t)GetInteger(index); }
int64_t SettingSFile::GetInteger64(unsigned int index) const { return GetInteger(index); }
uint64_t SettingSFile::GetUnsignedInteger64(unsigned int index) const { return (uint64_t)GetInteger(index); }
slibreal_t SettingSFile::GetScalar(unsigned int index) const {
    switch (this->dtype) {
        case DataType::INT32:  return (slibreal_t)this->data.i32->at(index);
        case DataType::INT64:  return (slibreal_t)this->data.i64->at(index);
        case DataType::UINT32: return (slibreal_t)this->data.u32->at(index);
        case DataType::UINT64: return (slibreal_t)this->data.u64->at(index);
        case DataType::REAL:   return this->data.r->at(index);

        default:
            throw ConfigurationException("The data in setting '%s' is not of numeric type.", this->name.c_str());
    }
}
string SettingSFile::GetString(unsigned int index) const {
    if (this->dtype == DataType::STRING)
        return this->data.s->at(index);
    else
        throw ConfigurationException("The data in setting '%s' is not of string type.", this->name.c_str());
}
vector<slibreal_t> SettingSFile::GetNumericVector() const {
    vector<slibreal_t> v(GetNumberOfValues());

    switch (this->dtype) {
        case DataType::INT32:
            for (int32_t a : *(this->data.i32))
                v.push_back((slibreal_t)a);
            break;
        case DataType::INT64:
            for (int64_t a : *(this->data.i64))
                v.push_back((slibreal_t)a);
            break;
        case DataType::UINT32:
            for (uint32_t a : *(this->data.u32))
                v.push_back((slibreal_t)a);
            break;
        case DataType::UINT64:
            for (uint64_t a : *(this->data.u64))
                v.push_back((slibreal_t)a);
            break;
        case DataType::REAL:
            for (slibreal_t a : *(this->data.r))
                v.push_back(a);
            break;

        default:
            throw ConfigurationException("The data in setting '%s' is of unrecognized type.", this->name.c_str());
    }

    return v;
}
const vector<string> SettingSFile::GetTextVector() const {
    if (this->dtype == DataType::STRING)
        return *(this->data.s);
    else
        throw ConfigurationException("The data in setting '%s' is not of string type.", this->name.c_str());
}

size_t SettingSFile::GetNumberOfValues() const {
    switch (this->dtype) {
        case DataType::BOOL:   return this->data.b->size();
        case DataType::INT32:  return this->data.i32->size();
        case DataType::INT64:  return this->data.i64->size();
        case DataType::UINT32: return this->data.u32->size();
        case DataType::UINT64: return this->data.u64->size();
        case DataType::REAL:   return this->data.r->size();
        case DataType::STRING: return this->data.s->size();

        default:
            throw ConfigurationException("The data in setting '%s' is of unrecognized type.", this->name.c_str());
    }
}

/************************
 * IS's                 *
 ************************/
bool SettingSFile::IsBool() const { return (this->dtype == DataType::BOOL); }
bool SettingSFile::IsBool(unsigned int) const { return IsBool(); }
bool SettingSFile::IsInteger() const {
    switch (this->dtype) {
        case DataType::INT32: case DataType::INT64:
        case DataType::UINT32: case DataType::UINT64:
        case DataType::REAL:
        return true;
        
        default: return false;
    }
}
bool SettingSFile::IsInteger(unsigned int) const { return IsInteger(); }
bool SettingSFile::IsInteger32() const {
    if (IsInteger()) {
        int64_t i = GetInteger();
        if (i <= std::numeric_limits<int32_t>::max() &&
            i >= std::numeric_limits<int32_t>::min())
            return true;
        else
            return false;
    } else
        return false;
}
bool SettingSFile::IsInteger32(unsigned int) const { return IsInteger32(); }
bool SettingSFile::IsInteger64() const { return (this->dtype == DataType::INT64); }
bool SettingSFile::IsInteger64(unsigned int) const { return IsInteger64(); }
bool SettingSFile::IsUnsignedInteger32() const { return IsUnsignedInteger32(0); }
bool SettingSFile::IsUnsignedInteger32(unsigned int index) const {
    // Make sure data is >= 0
    if (!IsUnsignedInteger64(index))
        return false;

    int64_t i = GetInteger(index);
    if (i <= std::numeric_limits<uint32_t>::max())
        return true;
    else
        return false;
}
bool SettingSFile::IsUnsignedInteger64() const { return IsUnsignedInteger64(0); }
bool SettingSFile::IsUnsignedInteger64(unsigned int index) const { 
    // Make sure data is >= 0
    switch (this->dtype) {
        case DataType::INT32:
            return (this->data.i32->at(index) >= 0);
        case DataType::INT64:
            return (this->data.i64->at(index) >= 0);
        case DataType::REAL:
            return (this->data.r->at(index) >= 0);

        case DataType::UINT32:
        case DataType::UINT64:
            return true;

        default:
            return false;
    }
}
bool SettingSFile::IsScalar() const {
    return (IsInteger() && this->GetNumberOfValues() == 1);
}
bool SettingSFile::IsScalar(unsigned int) const { return IsScalar(); }
bool SettingSFile::IsNumericVector() const { return IsInteger(); }
bool SettingSFile::IsNumericVector(unsigned int) const { return IsNumericVector(); }

/************************
 * SETTERS              *
 ************************/
void SettingSFile::SetValueBool(vector<bool> *b) {
    if (this->dtype != DataType::UNSPECIFIED)
        DeleteCurrent();

    this->dtype = DataType::BOOL;
    this->data.b = b;
}
void SettingSFile::SetValueInt32(vector<int32_t> *i) {
    if (this->dtype != DataType::UNSPECIFIED)
        DeleteCurrent();

    this->dtype = DataType::INT32;
    this->data.i32 = i;
}
void SettingSFile::SetValueInt64(vector<int64_t> *i) {
    if (this->dtype != DataType::UNSPECIFIED)
        DeleteCurrent();

    this->dtype = DataType::INT64;
    this->data.i64 = i;
}
void SettingSFile::SetValueUInt32(vector<uint32_t> *u) {
    if (this->dtype != DataType::UNSPECIFIED)
        DeleteCurrent();

    this->dtype = DataType::UINT32;
    this->data.u32 = u;
}
void SettingSFile::SetValueUInt64(vector<uint64_t> *u) {
    if (this->dtype != DataType::UNSPECIFIED)
        DeleteCurrent();

    this->dtype = DataType::UINT64;
    this->data.u64 = u;
}
void SettingSFile::SetValueReal(vector<slibreal_t> *s) {
    if (this->dtype != DataType::UNSPECIFIED)
        DeleteCurrent();

    this->dtype = DataType::REAL;
    this->data.r = s;
}
void SettingSFile::SetValueString(vector<string> *s) {
    if (this->dtype != DataType::UNSPECIFIED)
        DeleteCurrent();

    this->dtype = DataType::STRING;
    this->data.s = s;
}

