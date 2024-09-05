#include "fuzzyPID.h"

/*struct fuzzy {
    unsigned int input_num;// 输入变量的数量
    unsigned int output_num;//输出变量的数量
    unsigned int fo_type;//模糊系统的类型（通常表示模糊操作的类型）
    unsigned int *mf_type;//指向隶属度函数类型数组
    int *mf_params;// 指向模糊集参数的数组
    unsigned int df_type;// 指向模糊集参数的数组
    int *rule_base;//指向模糊规则库的数组
    float *output;// 指向输出结果的数组
};*/

struct fuzzy *fuzzy_init(unsigned int input_num, unsigned int output_num) {
    struct fuzzy *fuzzy_struct = (struct fuzzy *) malloc(sizeof(struct fuzzy));
    fuzzy_struct->input_num = input_num;
    fuzzy_struct->output_num = output_num;
    fuzzy_struct->mf_type = (unsigned int *) malloc((input_num + output_num) * sizeof(unsigned int));
#ifdef fuzzy_pid_rule_base_deep_copy
    fuzzy_struct->mf_params = (int *) malloc(4 * qf_default * sizeof(int));
    fuzzy_struct->rule_base = (int *) malloc(output_num * qf_default * qf_default * sizeof(int));
#endif
    fuzzy_struct->output = (float *) malloc(output_num * sizeof(float));
    return fuzzy_struct;
}

void delete_fuzzy(struct fuzzy *fuzzy_struct) {
    free(fuzzy_struct->mf_type);
    free(fuzzy_struct->output);
    free(fuzzy_struct);
}

void fuzzy_params_init(struct fuzzy *fuzzy_struct, unsigned int mf_type, unsigned int fo_type, unsigned int df_type,
                       int mf_params[], int rule_base[][qf_default])//初始化了一个 fuzzy 结构体中的一些参数。
{
    for (unsigned int i = 0; i < fuzzy_struct->input_num + fuzzy_struct->output_num; ++i) //循环代码用于初始化fuzzy_struct->mf_type数组。fuzzy_struct->input_num 和 fuzzy_struct->output_num 分别表示输入和输出的数量。
    {
        fuzzy_struct->mf_type[i] = mf_type;//这个循环将 mf_type 值赋给 mf_type 数组中的每一个元素。也就是说，所有输入和输出的隶属函数类型都被设置为 mf_type。
    }

    for (unsigned int i = 0; i < fuzzy_struct->output_num; ++i) //初始化 fuzzy_struct->output 数组。
    {
        fuzzy_struct->output[i] = 0;//将所有输出值初始化为 0。
    }

#ifdef fuzzy_pid_rule_base_deep_copy
    for (unsigned int j = 0; j < 4 * qf_default; ++j) {
         fuzzy_struct->mf_params[j] = mf_params[j];
     }

     for (unsigned int k = 0; k < fuzzy_struct->output_num * qf_default; ++k) {
         for (unsigned int i = 0; i < qf_default; ++i) {
             fuzzy_struct->rule_base[k * 7 + i] = rule_base[k][i];
         }
     }
#else
    fuzzy_struct->mf_params = mf_params;
    fuzzy_struct->rule_base = (int *) rule_base;
#endif

    fuzzy_struct->fo_type = fo_type;
    fuzzy_struct->df_type = df_type;
}

#define inverse(parameter) 1.0f/(float)parameter

// Gaussian membership function 高斯隶属度函数
float gaussmf(float x, float sigma, float c) {
    return expf(-powf(((x - c) / sigma), 2.0f));
}

// Generalized bell-shaped membership function 钟型隶属度函数
float gbellmf(float x, float a, float b, float c) {
    return inverse(1.0f + powf(fabsf((x - c) / a), 2.0f * b));
}

// Sigmoidal membership function s型隶属度函数
float sigmf(float x, float a, float c) {
    return inverse(1.0f + expf(a * (c - x)));
}

// Trapezoidal membership function T型隶属度函数
float trapmf(float x, float a, float b, float c, float d) {
    if (x >= a && x < b)
        return (x - a) / (b - a);
    else if (x >= b && x < c)
        return 1.0f;
    else if (x >= c && x <= d)
        return (d - x) / (d - c);
    else return 0.0f;
}

// Triangular membership function 三角隶属度函数
float trimf(float x, float a, float b, float c) {
    return trapmf(x, a, b, b, c);
}

// Z-shaped membership function Z型隶属度函数
float zmf(float x, float a, float b) {
    if (x <= a)
        return 1.0f;
    else if (x >= a && x <= (a + b) / 2.0f)
        return 1.0f - 2.0f * powf((x - a) / (b - a), 2.0f);
    else if (x >= (a + b) / 2.0f && x < b)
        return 2.0f * powf((x - b) / (b - a), 2.0f);
    else return 0;
}

// Membership function 隶属度函数
float mf(float x, unsigned int mf_type, int *params) //x是输入值，mf_type是隶属度函数类型，params：一个整型数组，包含用于计算隶属度函数的参数。具体的参数数量和含义取决于 mf_type 的值。
{
    switch (mf_type) //如果mf_type是
    {
        case 0://如果mf_type是0
            return gaussmf(x, params[0], params[1]);//返回高斯隶属度的计算值。模糊集合的元素、标准差 sigma、均值c。
        case 1:如果mf_type是1
            return gbellmf(x, params[0], params[1], params[2]);//返回钟型隶属度函数的计算值。模糊集合的元素、是宽度参数 a、是中心参数 b。
        case 2:如果mf_type是2
            return sigmf(x, params[0], params[2]);//返回s型隶属度函数的计算值。模糊集合的元素、是斜率参数 a、是中心参数 c。
        case 3:如果mf_type是3
            return trapmf(x, params[0], params[1], params[2], params[3]);//返回T型隶属度函数的计算值。模糊集合的元素、是梯形的左下角、 是梯形的左上角、是梯形的右上角、 是梯形的右下角
        case 4：如果mf_type是4
            return trimf(x,params[0], params[1], params[2]);
        case 5:如果mf_type是5
            return zmf(x, params[0], params[1]);//返回Z型隶属度函数的计算值。模糊集合的元素、是左侧的界限值、是右侧的界限值。
        default: // set triangular as default membership function将三角形设置为默认隶属函数
            return trimf(x, params[0], params[1], params[2]);////返回三角隶属度函数。模糊集合的元素、左端点、中间点（顶点）、和右端点。
    }
}

// Union operator
float or(float a, float b, unsigned int type)//计算并，ab是两个隶属度
{
    if (type == 1)//如果type == 1
    { // algebraic sum使用代数和的方法计算并集
        return a + b - a * b;
    } 
    else if (type == 2)//如果type == 2
    { // bounded sum使用有界和的方法计算并集
        return fminf(1, a + b);//结果限制在1之内，确保并集的最大值为1。
    } 
    else 
    { // fuzzy union模糊并集方法计算模糊并集
        return fmaxf(a, b);//取两者的最大值。
    }
}

// Intersection operator
float and(float a, float b, unsigned int type)//计算交。ab是两个隶属度
{
    if (type == 1) 如果type == 1
    { // algebraic product使用代数积方法计算模糊交集。 
        return a * b;//计算了两者的交集。
    }
    else if (type == 2)如果type == 2
    { // bounded product使用有界积方法计算模糊交集。确保结果不低于0，避免负值。
        return fmaxf(0, a + b - 1);//确保结果不低于0，避免负值。
    } 
    else 
    { // fuzzy intersection使用模糊交集方法计算模糊交集。
        return fminf(a, b);取两者的最小值。
    }
}

// Equilibrium operator
float equilibrium(float a, float b, float params)//计算平衡算子
{
    return powf(a * b, 1 - params) * powf(1 - (1 - a) * (1 - b), params);//powf(a * b, 1 - params)：考虑到交集的部分。powf(1 - (1 - a) * (1 - b), params)：考虑到并集的部分。
}

// Fuzzy operator
float fo(float a, float b, unsigned int type) {
    if (type < 3) {
        return and(a, b, type);
    } else if (type < 6) {
        return or(a, b, type - 3);
    } else {
        return equilibrium(a, b, 0.5f);
    }
}

// Mean of centers defuzzifier, only for two input multiple index
void moc(const float *joint_membership, const unsigned int *index, const unsigned int *count, struct fuzzy *fuzzy_struct) {

    float denominator_count = 0;
    float numerator_count[fuzzy_struct->output_num];
    for (unsigned int l = 0; l < fuzzy_struct->output_num; ++l) {
        numerator_count[l] = 0;
    }

    for (int i = 0; i < count[0]; ++i) {
        for (int j = 0; j < count[1]; ++j) {
            denominator_count += joint_membership[i * count[1] + j];
        }
    }

    for (unsigned int k = 0; k < fuzzy_struct->output_num; ++k) {
        for (unsigned int i = 0; i < count[0]; ++i) {
            for (unsigned int j = 0; j < count[1]; ++j) {
                numerator_count[k] += joint_membership[i * count[1] + j] *
                                      fuzzy_struct->rule_base[k * qf_default * qf_default + index[i] * qf_default +
                                                              index[count[0] + j]];
            }
        }
    }

#ifdef fuzzy_pid_debug_print
    printf("output:\n");
#endif
    for (unsigned int l = 0; l < fuzzy_struct->output_num; ++l) {
        fuzzy_struct->output[l] = numerator_count[l] / denominator_count;
#ifdef fuzzy_pid_debug_print
        printf("%f,%f,%f\n", numerator_count[l], denominator_count, fuzzy_struct->index[l]);
#endif
    }
}

// Defuzzifier
void df(const float *joint_membership, const unsigned int *output, const unsigned int *count, struct fuzzy *fuzzy_struct,
        int df_type) {
    if(df_type == 0)
        moc(joint_membership, output, count, fuzzy_struct);
    else {
        printf("Waring: No such of defuzzifier!\n");
        moc(joint_membership, output, count, fuzzy_struct);
    }
}

void fuzzy_control(float e, float de, struct fuzzy *fuzzy_struct) {
    float membership[qf_default * 2]; // Store membership
    unsigned int index[qf_default * 2]; // Store the index of each membership
    unsigned int count[2] = {0, 0};

    {
        int j = 0;
        for (int i = 0; i < qf_default; ++i) {
            float temp = mf(e, fuzzy_struct->mf_type[0], fuzzy_struct->mf_params + 4 * i);
            if (temp > 1e-4) {
                membership[j] = temp;
                index[j++] = i;
            }
        }

        count[0] = j;

        for (int i = 0; i < qf_default; ++i) {
            float temp = mf(de, fuzzy_struct->mf_type[1], fuzzy_struct->mf_params + 4 * i);
            if (temp > 1e-4) {
                membership[j] = temp;
                index[j++] = i;
            }
        }

        count[1] = j - count[0];
    }

#ifdef fuzzy_pid_debug_print
    printf("membership:\n");
    for (unsigned int k = 0; k < j; ++k) {
        printf("%f\n", membership[k]);
    }

    printf("index:\n");
    for (unsigned int k = 0; k < j; ++k) {
        printf("%d\n", index[k]);
    }

    printf("count:\n");
    for (unsigned int k = 0; k < 2; ++k) {
        printf("%d\n", count[k]);
    }
#endif

    if (count[0] == 0 || count[1] == 0) {
        for (unsigned int l = 0; l < fuzzy_struct->output_num; ++l) {
            fuzzy_struct->output[l] = 0;
        }
        return;
    }

    // Joint membership
    float joint_membership[count[0] * count[1]];

    for (int i = 0; i < count[0]; ++i) {
        for (int j = 0; j < count[1]; ++j) {
            joint_membership[i * count[1] + j] = fo(membership[i], membership[count[0] + j], fuzzy_struct->fo_type);
        }
    }

    df(joint_membership, index, count, fuzzy_struct, 0);
}

struct PID *raw_fuzzy_pid_init(float kp, float ki, float kd, float integral_limit, float dead_zone,
                               float feed_forward, float error_max, float delta_error_max, float delta_kp_max,
                               float delta_ki_max, float delta_kd_max, unsigned int mf_type, unsigned int fo_type,
                               unsigned int df_type, int mf_params[], int rule_base[][qf_default],
                               int output_min_value, int output_middle_value, int output_max_value) {
    struct PID *pid = (struct PID *) malloc(sizeof(struct PID));
    pid->kp = kp;
    pid->ki = ki;
    pid->kd = kd;

    pid->delta_kp_max = delta_kp_max;
    pid->delta_ki_max = delta_ki_max;
    pid->delta_kd_max = delta_kd_max;

    pid->delta_kp = 0;
    pid->delta_ki = 0;
    pid->delta_kd = 0;

    pid->error_max = error_max;
    pid->delta_error_max = delta_error_max;

    int output_count = 1;
    if (ki > 1e-4) {
        output_count += 1;
        if (kd > 1e-4)
            output_count += 1;
    }

    pid->fuzzy_struct = fuzzy_init(2, output_count);
    fuzzy_params_init(pid->fuzzy_struct, mf_type, fo_type, df_type, mf_params, rule_base);

    pid->last_error = 0;
    pid->current_error = 0;

    pid->intergral = 0;
    pid->intergral_limit = integral_limit;

    pid->dead_zone = dead_zone;
    pid->feed_forward = feed_forward;

    pid->output_max_value = output_max_value;
    pid->output_middle_value = output_middle_value;
    pid->output_min_value = output_min_value;

    return pid;
}

struct PID *fuzzy_pid_init(float *params, float delta_k, unsigned int mf_type, unsigned int fo_type,
                           unsigned int df_type, int mf_params[], int rule_base[][qf_default]) {
    return raw_fuzzy_pid_init(params[0], params[1], params[2], params[3], params[4], params[5], max_error,
                              max_delta_error, params[0] / delta_k, params[1] / delta_k, params[2] / delta_k, mf_type,
                              fo_type, df_type, mf_params,
                              rule_base, min_pwm_output, middle_pwm_output, max_pwm_output);
}

struct PID *raw_pid_init(float kp, float ki, float kd, float integral_limit, float dead_zone,
                         float feed_forward, float linear_adaptive_kp, float error_max, float delta_error_max,
                         int output_min_value, int output_middle_value, int output_max_value) {
    struct PID *pid = (struct PID *) malloc(sizeof(struct PID));
    pid->kp = kp;
    pid->ki = ki;
    pid->kd = kd;

    pid->delta_kp_max = 0;
    pid->delta_ki_max = 0;
    pid->delta_kd_max = 0;

    pid->delta_kp = 0;
    pid->delta_ki = 0;
    pid->delta_kd = 0;

    pid->error_max = error_max;
    pid->delta_error_max = delta_error_max;

    pid->fuzzy_struct = NULL;

    pid->last_error = 0;
    pid->current_error = 0;

    pid->intergral = 0;
    pid->intergral_limit = integral_limit;

    pid->dead_zone = dead_zone;
    pid->feed_forward = feed_forward;

    pid->output_max_value = output_max_value;
    pid->output_middle_value = output_middle_value;
    pid->output_min_value = output_min_value;

    pid->linear_adaptive_kp = linear_adaptive_kp;
    return pid;
}

struct PID *pid_init(float *params) {
    return raw_pid_init(params[0], params[1], params[2], params[3], params[4], params[5], params[6], max_error,
                        max_delta_error, min_pwm_output, middle_pwm_output, max_pwm_output);
}

int round_user(float parameter) {
    if ((int) (parameter * 10.0) % 10 >= 5)
        return parameter + 1;
    else
        return parameter;
}

int limit(int value, int max_limit, int min_limit) {
    if (value > max_limit)
        return max_limit;
    if (value < min_limit)
        return min_limit;
    return value;
}

float fuzzy_pid_control(float real, float idea, struct PID *pid) {
    pid->last_error = pid->current_error;
    pid->current_error = idea - real;
    float delta_error = pid->current_error - pid->last_error;
#ifdef fuzzy_pid_dead_zone
    if (pid->current_error < pid->dead_zone && pid->current_error > -pid->dead_zone) {
        pid->current_error = 0;
    } else {
        if (pid->current_error > pid->dead_zone)
            pid->current_error = pid->current_error - pid->dead_zone;
        else {
            if (pid->current_error < -pid->dead_zone)
                pid->current_error = pid->current_error + pid->dead_zone;
        }
    }
#endif
    fuzzy_control(pid->current_error / pid->error_max * 3.0f, delta_error / pid->delta_error_max * 3.0f,
                  pid->fuzzy_struct);

    pid->delta_kp = pid->fuzzy_struct->output[0] / 3.0f * pid->delta_kp_max + pid->kp;

    if (pid->fuzzy_struct->output_num >= 2)
        pid->delta_ki = pid->fuzzy_struct->output[1] / 3.0f * pid->delta_ki_max;
    else pid->delta_ki = 0;

    if (pid->fuzzy_struct->output_num >= 3)
        pid->delta_kd = pid->fuzzy_struct->output[2] / 3.0f * pid->delta_kd_max;
    else pid->delta_ki = 0;

#ifdef fuzzy_pid_debug_print
    printf("kp : %f, ki : %f, kd : %f\n", kp, ki, kd);
#endif

    pid->intergral += (pid->ki + pid->delta_ki) * pid->current_error;
#ifdef fuzzy_pid_integral_limit
    if (pid->intergral > pid->intergral_limit)
        pid->intergral = pid->intergral_limit;
    else {
        if (pid->intergral < -pid->intergral_limit)
            pid->intergral = -pid->intergral_limit;
    }
#endif
    pid->output = (pid->kp + pid->delta_kp) * pid->current_error + pid->intergral +
                  (pid->kd + pid->delta_kd) * (pid->current_error - pid->last_error);
    pid->output += pid->feed_forward * (float) idea;
    return pid->output;
}


float pid_control(float real, float idea, struct PID *pid) {
    pid->last_error = pid->current_error;
    pid->current_error = idea - real;

#ifdef pid_dead_zone
    if (pid->current_error < pid->dead_zone && pid->current_error > -pid->dead_zone) {
        pid->current_error = 0;
    } else {
        if (pid->current_error > pid->dead_zone)
            pid->current_error = pid->current_error - pid->dead_zone;
        else {
            if (pid->current_error < -pid->dead_zone)
                pid->current_error = pid->current_error + pid->dead_zone;
        }
    }
#endif

#ifdef pid_debug_print
    printf("kp : %f, ki : %f, kd : %f\n", kp, ki, kd);
#endif

    pid->intergral += (pid->ki) * pid->current_error;
#ifdef pid_integral_limit
    if (pid->intergral > pid->intergral_limit)
        pid->intergral = pid->intergral_limit;
    else {
        if (pid->intergral < -pid->intergral_limit)
            pid->intergral = -pid->intergral_limit;
    }
#endif

    float linear_adaptive_kp = 1;
    if (pid->linear_adaptive_kp > 1e-4)
        linear_adaptive_kp =
                (1 - pid->linear_adaptive_kp) * pid->current_error / pid->error_max + pid->linear_adaptive_kp;

    pid->output = pid->kp * linear_adaptive_kp * pid->current_error + pid->intergral +
                  (pid->kd) * (pid->current_error - pid->last_error);
    pid->output += pid->feed_forward * (float) idea;
    return pid->output;
}


void delete_pid(struct PID *pid) {
    if (pid->fuzzy_struct != NULL) {
        delete_fuzzy(pid->fuzzy_struct);
    }
    free(pid);
}

void delete_pid_vector(struct PID **pid_vector, unsigned int count) {
    for (unsigned int i = 0; i < count; ++i) {
        delete_pid(pid_vector[i]);
    }
    free(pid_vector);
}

struct PID **pid_vector_init(float params[][pid_params_count], unsigned int count) {
    struct PID **pid = (struct PID **) malloc(sizeof(struct PID *) * count);
    for (unsigned int i = 0; i < count; ++i) {
        pid[i] = pid_init(params[i]);
    }
    return pid;
}

struct PID **
fuzzy_pid_vector_init(float params[][pid_params_count], float delta_k, unsigned int mf_type, unsigned int fo_type,
                      unsigned int df_type, int *mf_params, int rule_base[][qf_default],
                      unsigned int count) {
    struct PID **pid = (struct PID **) malloc(sizeof(struct PID *) * count);
    for (unsigned int i = 0; i < count; ++i) {
        pid[i] = fuzzy_pid_init(params[i], delta_k, mf_type, fo_type, df_type, mf_params, rule_base);
    }
    return pid;
}

int direct_control(int zero_value, int offset_value, bool direct) {
    if (direct == true) {
        return zero_value + offset_value;
    } else {
        return zero_value - offset_value;
    }
}

int fuzzy_pid_motor_pwd_output(float real, float idea, bool direct, struct PID *pid) {
    return limit(direct_control(pid->output_middle_value, fuzzy_pid_control(real, idea, pid), direct),
                 pid->output_max_value, pid->output_min_value);
}

int pid_motor_pwd_output(float real, float idea, bool direct, struct PID *pid) {
    return limit(direct_control(pid->output_middle_value, pid_control(real, idea, pid), direct), pid->output_max_value,
                 pid->output_min_value);
}

