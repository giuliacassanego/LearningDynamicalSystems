# AsLorenzoSaid - Learning Dynamical Systems

This project focuses on **Learning Dynamical Systems**, specifically for **Optical Flow Estimation** using filtered data (EKF, UKF). The implementation is based on the [PX4-Autopilot v1.16.0](https://github.com/PX4/PX4-Autopilot).

## Project Structure

- **`Data/`**: Contains optical flow datasets. Note that raw and processed data are excluded from version control via `.gitignore`.
  - `estimator_optical_flow`: Velocity in NED frame using EKF.
  - `vehicle_optical_flow`: Processed optical flow data (aligned using sampling time).
  - `sensor_optical_flow`: Raw sensor data.
- **`filters/`**: Core implementations of various Kalman Filters:
  - `EKF.m` / `REKF.m`: Extended Kalman Filter / Robust EKF.
  - `UKF.m` / `RUKF.m`: Unscented Kalman Filter / Robust UKF.
  - `func_f.m`, `func_h.m`: System dynamics and measurement functions.
- **`project/`**: Contains specific analysis and verification tasks.
- **`presentation/`**: Presentation materials for the course.

## Project Tasks

The project is organized into several tasks located in the `project/` directory:

- **Task 1**: `optical_flow_model.m` - Implementation of the optical flow state-space model.
- **Task 2**: `task2.m` - Accuracy verification of the optical flow model.

## Data Handling

To avoid uploading large datasets and logs:
- `.gitignore` is used to exclude `.mat` files and `log_*/` directories within `Data/mat/`.
- Refer to `Data/README.txt` for more details on the data structure.

## Requirements
- MATLAB
- (Optional) Navigation Toolbox for conversions (e.g., `quat2rotm`).