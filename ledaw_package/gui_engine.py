from PySide6.QtWidgets import QMainWindow, QApplication, QWidget, QFrame, QScrollArea, QDialog, QTextEdit, QProgressBar, QDialogButtonBox, QFileDialog, QVBoxLayout, QLabel, QPushButton, QLineEdit, QTabWidget, QFormLayout, QComboBox, QCheckBox, QMessageBox, QHBoxLayout, QSizePolicy, QSpacerItem
from PySide6.QtGui import QIcon
from PySide6.QtCore import Qt, QTimer, Signal, Slot
import matplotlib.pyplot as plt
from .job_engine import JobWorker, PlotJobWorker
import shutil
import sys
import os


if sys.platform.lower().startswith('win'):
    operating_system = 'windows'
elif sys.platform == 'darwin':
    operating_system = 'mac'
else:
    operating_system = 'linux'


class LEDAWApp(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle('LEDAW - LED Analysis Wizard')

        # Set icon path based on whether the app is bundled with PyInstaller or running in a normal Python environment
        if getattr(sys, 'frozen', False):  # If the application is bundled by PyInstaller
            base_path = sys._MEIPASS  # This is where PyInstaller stores temporary files
        else:
            base_path = os.path.dirname(os.path.abspath(sys.argv[0]))  # Normal environment (not bundled)

        # Set window icon based on the operating system
        try:
            if sys.platform.startswith('win'):
                icon_path = os.path.join(base_path, 'ledaw_package', 'favicon.ico')
            elif sys.platform == 'darwin':  # macOS
                icon_path = os.path.join(base_path, 'ledaw_package', 'favicon.icns')
            else:  # Assume Linux or other OS
                icon_path = os.path.join(base_path, 'ledaw_package', 'favicon.png')

            # Set the window icon
            self.setWindowIcon(QIcon(icon_path))

        except Exception as e:
            print(f"Could not set the window icon: {e}")
            # Continue without crashing if the icon cannot be set

        self.setGeometry(100, 100, 600, 400)

        # Initialize flags for the tabs
        self.nbody_tab = None
        self.twobody_tab = None
        self.plot_tab = None
        self.created_tabs = [] 
        self.perform_nbody = False
        self.perform_twobody = False
        self.perform_plot = False
        self.perform_cps = False
        self.perform_cbs = False
        self.bypass_validation = False
        self.nbody_locked = False
        self.twobody_locked = False
        self.plot_locked = False
        self.show_diag_cells_for_fp_led = False
        self.delete_old_plots = False
    
        # Create Tabs
        self.tabs = QTabWidget()
        self.create_home_tab()  # Always create Home tab first
        self.create_manage_session_tab()  # Ensure this tab is always present
        self.create_about_tab()  # Ensure this tab is always present

        # Main Layout
        layout = QVBoxLayout()
        layout.addWidget(self.tabs)

        # Create a container widget to hold the layout
        self.container = QWidget(self)
        self.container.setLayout(layout)

        # Set this container widget as the central widget for the main window
        self.setCentralWidget(self.container)

        # Initialize variables
        self.ledaw_out_root = None
        self.is_manual_directory = False  # Flag to track if ledaw output directory was manually set
        self.conversion_factor = 627.5095  # Default to kcal/mol
        self.relabel_mapping = {}
        self.job_name = None
        self.f_ref = None
        self.f_corr = None
        self.F = None
        self.last_selected_dir = './'
        self.plot_unlock_count = 0  # Initialize a counter for Plot tab unlocks

        # Connect the 'editingFinished' event to trigger the update after typing is done
        self.output_dir_input.editingFinished.connect(self.update_output_directory)


    @Slot(str, str)
    def show_error_message(self, title, message):
        QMessageBox.critical(self, title, message)


    def normalize_path_gui(self, path):
        """Normalize the path by converting backslashes to forward slashes and removing trailing slashes."""
        # Replace backslashes with forward slashes for consistency
        normalized_path = path.replace("\\", "/")
        return normalized_path.rstrip("/")  # remove any trailing slashes


    def update_filename(self, file_list, new_value, index):
        """Ensure the list size and update the filename with normalized path."""
        # Ensure the list size is adequate
        self.ensure_list_size(file_list, index)

        # Normalize the path and update the list
        file_list[index] = self.normalize_path_gui(new_value.strip() or '')


    def change_output_directory(self):
        """Allow the user to select a different output directory."""
        selected_dir = QFileDialog.getExistingDirectory(self, "Select Output Directory", "./LEDAW-OUT")
        if selected_dir:
            self.ledaw_out_root = self.normalize_path_gui(selected_dir)
            self.is_manual_directory = True  # Mark as manually updated
            # Update the text box with the selected directory
            self.output_dir_input.blockSignals(True)
            self.output_dir_input.setText(self.ledaw_out_root)
            self.output_dir_input.blockSignals(False)


    def update_output_directory(self):
        """Update the LEDAW output directory when the user manually edits the directory path."""
        # Normalize the manually entered directory path
        self.ledaw_out_root = self.normalize_path_gui(self.output_dir_input.text())
        self.is_manual_directory = True  # Mark as manually updated


    def showEvent(self, event):
        # Ensure the icon is properly set once the window is shown
        try:
            if getattr(sys, 'frozen', False):  # If running as a PyInstaller bundle
                base_path = sys._MEIPASS
            else:
                base_path = os.path.dirname(os.path.abspath(sys.argv[0]))

            # Set the icon based on the operating system
            if sys.platform.startswith('win'):
                icon_path = os.path.join(base_path, 'ledaw_package', 'favicon.ico')
            elif sys.platform == 'darwin':  # macOS
                icon_path = os.path.join(base_path, 'ledaw_package', 'favicon.icns')
            else:  # Assume Linux or other OS
                icon_path = os.path.join(base_path, 'ledaw_package', 'favicon.png')

            # Set the window icon
            self.setWindowIcon(QIcon(icon_path))

        except Exception as e:
            print(f"Could not set the window icon: {e}")
            # Continue without crashing if the icon cannot be set

        # Call the base class's showEvent method
        super().showEvent(event)


    def show_job_name_info(self):
        msg = QMessageBox()
        msg.setWindowTitle("Job Name Info")

        try:
            # If running under PyInstaller, use sys._MEIPASS for the path
            if getattr(sys, 'frozen', False):  # Running as a PyInstaller bundle
                base_path = sys._MEIPASS
            else:
                base_path = os.path.dirname(os.path.abspath(sys.argv[0]))

            # Set the icon based on the operating system
            if sys.platform.startswith('win'):
                icon_path = os.path.join(base_path, 'ledaw_package', 'favicon.ico')
            elif sys.platform == 'darwin':  # macOS
                icon_path = os.path.join(base_path, 'ledaw_package', 'favicon.icns')
            else:  # Assume Linux or other OS
                icon_path = os.path.join(base_path, 'ledaw_package', 'favicon.png')

            # Set the message box icon
            msg.setWindowIcon(QIcon(icon_path))

        except Exception as e:
            print(f"Could not set the message box icon: {e}")
            # Continue without crashing if the icon fails to load

        # Set the informative text for the message box
        msg.setText(
            "Example: WATER-HEXAMER\n\n"
            "Note: By default, a subdirectory with the specified job name will be created in\n"
            "./LEDAW-OUT/ to store the LEDAW results. If the specified subdirectory already exists, "
            "LEDAW will overwrite its contents, and existing files/subdirectories may interfere with the LEDAW run. "
            "It is strongly recommended to use a unique job name that does not match existing subdirectory names "
            "under ./LEDAW-OUT/ to avoid potential conflicts.\n\n"
            "If needed, the default ./LEDAW-OUT/ path can be changed by filling the next field."
        )

        # Show the message box
        msg.exec()


    def show_method_info(self):
        msg = QMessageBox()
        msg.setWindowTitle("Method Info")

        try:
            # If running under PyInstaller, use sys._MEIPASS for the path
            if getattr(sys, 'frozen', False):  # Running as a PyInstaller bundle
                base_path = sys._MEIPASS
            else:
                base_path = os.path.dirname(os.path.abspath(sys.argv[0]))

            # Set the icon based on the operating system
            if sys.platform.startswith('win'):
                icon_path = os.path.join(base_path, 'ledaw_package', 'favicon.ico')
            elif sys.platform == 'darwin':  # macOS
                icon_path = os.path.join(base_path, 'ledaw_package', 'favicon.icns')
            else:  # Assume Linux or other OS
                icon_path = os.path.join(base_path, 'ledaw_package', 'favicon.png')

            # Set the message box icon
            msg.setWindowIcon(QIcon(icon_path))

        except Exception as e:
            print(f"Could not set the message box icon: {e}")
            # Continue without crashing if the icon fails to load
    
        msg.setText("For several reasons, the method name is not directly extracted from ORCA output files. "
                    "LED terms will be collected from the output files according to the chosen method.")
        msg.exec()


    def create_home_tab(self):
        
        self.home_tab = QWidget() 
        self.layout = QFormLayout()
    
        # Job Name input
        self.job_name_input = QLineEdit()
        job_name_info_button = QPushButton('i')
        job_name_info_button.setFixedSize(25, 25)
        job_name_info_button.setStyleSheet("border-radius: 12px; background-color: lightgray; color: black; font-weight: bold; font-size: 13px;")
        job_name_info_button.clicked.connect(self.show_job_name_info)

        job_name_layout = QHBoxLayout()
        job_name_layout.addWidget(self.job_name_input)
        job_name_layout.addWidget(job_name_info_button)
        self.layout.addRow("Assign a name to LEDAW job:", job_name_layout)

        # Validate job name
        self.job_name_input.textChanged.connect(self.validate_job_name)

        # Set LEDAW output root directory
        self.output_dir_label = QLabel("LEDAW outputs will be written to:")
        self.output_dir_input = QLineEdit()
        self.output_dir_input.textChanged.connect(self.update_output_directory)
        self.change_dir_button = QPushButton("Change Directory")
        self.change_dir_button.clicked.connect(self.change_output_directory)

        # Add output directory layout
        output_dir_layout = QHBoxLayout()
        output_dir_layout.addWidget(self.output_dir_label)
        output_dir_layout.addWidget(self.output_dir_input)
        output_dir_layout.addWidget(self.change_dir_button)

        self.layout.addRow(output_dir_layout)

        # Method Selection (set DLPNO-CCSD(T) as default)
        self.method = QComboBox()
        self.method.addItems(["DLPNO-CCSD(T)", "DLPNO-CCSD", "HFLD"])
        self.method.setCurrentText("DLPNO-CCSD(T)")
        self.method.currentTextChanged.connect(self.toggle_rhf_checkbox)

        method_info_button = QPushButton('i')
        method_info_button.setFixedSize(25, 25)
        method_info_button.setStyleSheet("border-radius: 12px; background-color: lightgray; color: black; font-weight: bold; font-size: 13px;")
        method_info_button.clicked.connect(self.show_method_info)

        method_layout = QHBoxLayout()
        method_layout.addWidget(self.method)
        method_layout.addWidget(method_info_button)
        self.layout.addRow("Choose a method:", method_layout)

        # RHF checkbox (initially hidden)
        self.use_rhf_checkbox = QCheckBox("Use RHF energy as reference in HFLD instead of E(0) for subsystems?")
        self.use_rhf_checkbox.setVisible(False)

        # Create an RHF layout
        self.rhf_layout = QHBoxLayout()
        self.rhf_layout.addWidget(self.use_rhf_checkbox)

        self.layout.addRow(self.rhf_layout)

        # Interaction Energy Unit Selection
        self.unit_label = QLabel("In which unit do you want to calculate interaction energies?")
        self.unit_choice = QComboBox()
        self.unit_choice.addItems(["kcal/mol", "kJ/mol"])
        self.unit_choice.setCurrentText("kcal/mol")
        self.unit_choice.currentTextChanged.connect(self.set_conversion_factor)
        self.layout.addRow(self.unit_label, self.unit_choice)

        # Add questions for N-body and Two-body LED
        self.nbody_checkbox = QCheckBox("Do you want N-body LED?")
        nbody_info_button = QPushButton('i')
        nbody_info_button.setFixedSize(25, 25)  # Minimal size
        nbody_info_button.setStyleSheet("border-radius: 12px; background-color: lightgray; color: black; font-weight: bold; font-size: 13px;")
        nbody_info_button.clicked.connect(self.show_nbody_info)

        self.twobody_checkbox = QCheckBox("Do you want two-body LED?")
        twobody_info_button = QPushButton('i')
        twobody_info_button.setFixedSize(25, 25)  # Minimal size
        twobody_info_button.setStyleSheet("border-radius: 12px; background-color: lightgray; color: black; font-weight: bold; font-size: 13px;")
        twobody_info_button.clicked.connect(self.show_twobody_info)

        # Layout for N-body and Two-body LED
        row1_layout = QHBoxLayout()
        row1_layout.addWidget(self.nbody_checkbox)
        row1_layout.addWidget(nbody_info_button)
        self.layout.addRow(row1_layout)

        row2_layout = QHBoxLayout()
        row2_layout.addWidget(self.twobody_checkbox)
        row2_layout.addWidget(twobody_info_button)
        self.layout.addRow(row2_layout)

        # CPS Extrapolation Checkbox
        self.cps_checkbox = QCheckBox("Do you want CPS extrapolation?")
        self.cps_checkbox.toggled.connect(self.set_cps_f_value)  # Connect to set_cps_f_value function
        cps_info_button = QPushButton('i')
        cps_info_button.setFixedSize(25, 25)
        cps_info_button.setStyleSheet("border-radius: 12px; background-color: lightgray; color: black; font-weight: bold; font-size: 13px;")
        cps_info_button.clicked.connect(self.show_cps_info)

        cps_layout = QHBoxLayout()
        cps_layout.addWidget(self.cps_checkbox)
        cps_layout.addWidget(cps_info_button)
        self.layout.addRow(cps_layout)

        # CBS Extrapolation Checkbox and Method Selection
        self.cbs_checkbox = QCheckBox("Do you want CBS extrapolation?")
        cbs_info_button = QPushButton('i')
        cbs_info_button.setFixedSize(25, 25)
        cbs_info_button.setStyleSheet("border-radius: 12px; background-color: lightgray; color: black; font-weight: bold; font-size: 13px;")
        cbs_info_button.clicked.connect(self.show_cbs_info)
        self.cbs_checkbox.toggled.connect(self.toggle_cbs_options)

        cbs_layout = QHBoxLayout()
        cbs_layout.addWidget(self.cbs_checkbox)
        cbs_layout.addWidget(cbs_info_button)
        self.layout.addRow(cbs_layout)

        # CBS Method Selection (Label, ComboBox, and F coefficients)
        self.cbs_label = QLabel("Select CBS method:")
        self.cbs_option = QComboBox()
        self.cbs_option.addItems(["CBS(3/4)", "CBS(2/3)", "Other"])
        self.cbs_option.setFixedWidth(160)
        self.cbs_option.currentIndexChanged.connect(self.set_cbs_defaults)  # Initialize the default CBS values

        self.f_ref_input = QLineEdit()
        self.f_ref_input.setPlaceholderText("F coefficient for reference part")
        self.f_corr_input = QLineEdit()
        self.f_corr_input.setPlaceholderText("F coefficient for correlation part")

        # CBS validation error label (initially hidden)
        self.cbs_error_label = QLabel('')
        self.cbs_error_label.setStyleSheet("color: red;")
        self.cbs_error_label.setVisible(False)

        # Add CBS method, dropdown, and F coefficient fields all in one line
        cbs_method_layout = QHBoxLayout()
        cbs_method_layout.addWidget(self.cbs_label)
        cbs_method_layout.addWidget(self.cbs_option)
        cbs_method_layout.addWidget(self.f_ref_input)
        cbs_method_layout.addWidget(self.f_corr_input)

        # Add CBS error label below the F inputs
        self.layout.addRow(cbs_method_layout)
        self.layout.addRow(self.cbs_error_label)  # Error message row

        # Initially hide F inputs
        self.cbs_label.setVisible(False)
        self.cbs_option.setVisible(False)
        self.f_ref_input.setVisible(False)
        self.f_corr_input.setVisible(False)        
        
        # Ensure alignment from left
        cbs_method_layout.setAlignment(Qt.AlignmentFlag.AlignLeft)

        # Initially hide options to CBS method
        self.cbs_label.setVisible(False)
        self.cbs_option.setVisible(False)
        self.f_ref_input.setVisible(False)
        self.f_corr_input.setVisible(False)

        # LED Heat Map Plot
        self.plot_checkbox = QCheckBox("Do you want to plot LED heat maps?")
        plot_info_button = QPushButton('i')
        plot_info_button.setFixedSize(25, 25)  # Minimal size
        plot_info_button.setStyleSheet("border-radius: 12px; background-color: lightgray; color: black; font-weight: bold; font-size: 13px;")
        plot_info_button.clicked.connect(self.show_plot_info)

        row3_layout = QHBoxLayout()
        row3_layout.addWidget(self.plot_checkbox)
        row3_layout.addWidget(plot_info_button)
        self.layout.addRow(row3_layout)

        # Button to show parameters
        show_params_button = QPushButton("Confirm Parameters to Proceed")
        show_params_button.clicked.connect(self.show_parameters)
        show_params_button.setFixedHeight(35)
        show_params_button.setStyleSheet("QPushButton {background-color: #ebedef;} QPushButton:hover {background-color: paleturquoise;}")
        self.layout.addRow(show_params_button)

        # Connect checkboxes to toggle the flags
        self.nbody_checkbox.toggled.connect(lambda: setattr(self, 'perform_nbody', self.nbody_checkbox.isChecked()))
        self.twobody_checkbox.toggled.connect(lambda: setattr(self, 'perform_twobody', self.twobody_checkbox.isChecked()))
        self.plot_checkbox.toggled.connect(lambda: setattr(self, 'perform_plot', self.plot_checkbox.isChecked()))

        self.home_tab.setLayout(self.layout)
        self.tabs.addTab(self.home_tab, "Home")


    def toggle_rhf_checkbox(self):
        if self.method.currentText() == "HFLD":
            if not self.use_rhf_checkbox.isVisible():
                self.use_rhf_checkbox.setVisible(True)
                self.layout.insertRow(2, self.rhf_layout)  # Insert below the method selection
        else:
            self.use_rhf_checkbox.setVisible(False)


    def set_conversion_factor(self):
        """Convert and set the conversion factor to float."""
        if self.unit_choice.currentText() == "kcal/mol":
            self.conversion_factor = 627.5095
        elif self.unit_choice.currentText() == "kJ/mol":
            self.conversion_factor = 2625.5
        # Ensure conversion_factor is treated as float for future calculations
        self.conversion_factor = float(self.conversion_factor)


    def set_cps_f_value(self):
        """Convert and set CPS F value to float if selected."""
        if self.cps_checkbox.isChecked():
            self.F = 1.5
        else:
            self.F = None
        # Ensure F is treated as float for future calculations
        if self.F is not None:
            self.F = float(self.F)


    def toggle_directory_selection(self):
        """Show or hide the directory selection option based on Yes/No dropdown."""
        if self.change_dir_dropdown.currentText() == "Yes":
            self.change_dir_button.setVisible(True)
        else:
            self.change_dir_button.setVisible(False)


    def toggle_cbs_options(self):
        """Show/hide CBS method options and F coefficient inputs."""
        if self.cbs_checkbox.isChecked():
            # Show CBS method selection and default CBS values when checked
            self.cbs_label.setVisible(True)
            self.cbs_option.setVisible(True)
            self.set_cbs_defaults()  # Initialize default CBS values when checked
        else:
            # Hide CBS method selection and F coefficient inputs
            self.cbs_label.setVisible(False)
            self.cbs_option.setVisible(False)
            self.f_ref_input.setVisible(False)
            self.f_corr_input.setVisible(False)


    def show_scrollable_confirmation_dialog(self, title, content):
        """Display a confirmation dialog with a scrollable area for the content."""
        dialog = QDialog(self)
        dialog.setWindowTitle(title)

        # Create a vertical layout for the dialog
        layout = QVBoxLayout(dialog)

        # Create a scroll area and a widget inside it for the content
        scroll_area = QScrollArea(dialog)
        scroll_area.setWidgetResizable(True)  # Allow the scroll area to resize with the window

        # Create a widget to hold the content (this will be the scrollable area)
        scroll_content = QWidget()
        scroll_layout = QVBoxLayout(scroll_content)

        # Add the content (the overview of files) to the scrollable layout
        content_label = QLabel(content)
        content_label.setWordWrap(True)  # Enable word wrap for long content
        content_label.setAlignment(Qt.AlignmentFlag.AlignLeft)  # Align the content to the left
        scroll_layout.addWidget(content_label)

        # Set the scrollable content widget as the scroll area's widget
        scroll_area.setWidget(scroll_content)

        # Add the scroll area to the dialog's layout
        layout.addWidget(scroll_area)

        # Add confirmation buttons
        button_layout = QHBoxLayout()
        yes_button = QPushButton("Yes", dialog)
        no_button = QPushButton("No", dialog)

        yes_button.clicked.connect(lambda: dialog.accept())  # Close dialog with accept
        no_button.clicked.connect(lambda: dialog.reject())  # Close dialog with reject

        # Add buttons to layout
        button_layout.addWidget(yes_button)
        button_layout.addWidget(no_button)
        layout.addLayout(button_layout)

        # Show the dialog
        result = dialog.exec()

        return result == QDialog.DialogCode.Accepted


    def set_cbs_defaults(self):
        """Set default F values based on CBS method and show/hide F coefficient inputs."""
        cbs_method = self.cbs_option.currentText()
        if cbs_method == "CBS(3/4)":
            self.f_ref_input.setText("1.30130392358216")
            self.f_corr_input.setText("1.71188948355961")
            self.f_ref_input.setVisible(False)
            self.f_corr_input.setVisible(False)
        elif cbs_method == "CBS(2/3)":
            self.f_ref_input.setText("1.32521623414431")
            self.f_corr_input.setText("1.58433631967161")
            self.f_ref_input.setVisible(False)
            self.f_corr_input.setVisible(False)
        elif cbs_method == "Other":
            self.f_ref_input.clear()
            self.f_corr_input.clear()
            self.f_ref_input.setVisible(True)
            self.f_corr_input.setVisible(True)


    def validate_f_values_immediate(self):
        """Validate the F values immediately after the user inputs them."""
        try:
            # Attempt to convert the input to float
            f_ref = float(self.f_ref_input.text())
            f_corr = float(self.f_corr_input.text())

            # If both are valid, set them as floats and hide any error messages
            self.F_ref_cbs = f_ref
            self.F_corr_cbs = f_corr
            self.cbs_error_label.setVisible(False)  # Hide error label if validation is successful
            return True  # Return True indicating successful validation

        except ValueError:
            # Show a pop-up window with the error message
            QMessageBox.warning(self, "Invalid Input", "Please enter valid numerical values for F coefficients.")
            self.cbs_error_label.setVisible(True)
            return False  # Return False indicating validation failure


    def set_cps_f_value(self):
        """Set F to 1.5 if CPS is selected, otherwise set it to None."""
        if self.cps_checkbox.isChecked():
            self.F = float(1.5)
        else:
            self.F = None
        return self.F


    def validate_f_values(self):
        if not self.cbs_checkbox.isChecked():
            return None, None  # Do not validate if CBS is not selected
        try:
            self.f_ref = float(self.f_ref_input.text())
            self.f_corr = float(self.f_corr_input.text())
            return self.f_ref, self.f_corr
        except ValueError:
            QMessageBox.warning(self, "Invalid Input", "Please enter valid numerical values for F coefficients.")
            return None, None


    def validate_method(self):
        method = self.method.currentText()
        if method not in ["DLPNO-CCSD(T)", "DLPNO-CCSD", "HFLD"]:
            QMessageBox.warning(self, "Invalid Method", "Please select a valid method.")
            return None
        return method


    def validate_job_name(self):
        """Validates the job name and sets the output directory."""
        job_name = self.job_name_input.text().strip()

        if not job_name:  # If job name is empty, set to None
            self.job_name = None
            self.ledaw_out_root = None
            self.output_dir_input.setText('')
        elif " " in job_name:
            QMessageBox.warning(self, "Invalid Job Name", "Job name cannot contain spaces.")
            self.job_name = None
            self.ledaw_out_root = None
            self.output_dir_input.setText('')
        else:
            self.job_name = job_name
            self.ledaw_out_root = rf"./LEDAW-OUT/{self.job_name}"
            self.output_dir_input.setText(self.ledaw_out_root)


    def check_selections(self):
        """Ensure that at least one of N-body, Two-body, and Plot is selected."""
        if self.bypass_validation:
            return True  # Skip validation when bypass is active

        # Standard validation: ensure at least one of N-body, Two-body, or Plot is selected
        if not (self.nbody_checkbox.isChecked() or self.twobody_checkbox.isChecked() or self.plot_checkbox.isChecked()):
            QMessageBox.warning(self, "Invalid Selection", 
                                "At least one of N-body LED, Two-body LED, and Plot LED maps must be selected.")
            return False

        return True


    def get_rhf_selection(self):
        return self.use_rhf_checkbox.isChecked() if self.use_rhf_checkbox.isVisible() else None


    def format_f_values_for_display(self, f_ref, f_corr):
        """Format F values for display based on CBS method."""
        if self.cbs_option.currentText() in ["CBS(3/4)", "CBS(2/3)"]:
            # Fixed 6 decimal places for CBS(3/4) and CBS(2/3)
            f_ref_str = f"{float(f_ref):.6f}"
            f_corr_str = f"{float(f_corr):.6f}"
        else:
            # Dynamically adjust precision for 'Other'
            ref_decimals = len(f_ref.split('.')[-1]) if '.' in f_ref else 0
            corr_decimals = len(f_corr.split('.')[-1]) if '.' in f_corr else 0
            max_decimals = max(ref_decimals, corr_decimals)

            f_ref_str = f"{float(f_ref):.{max_decimals}f}"
            f_corr_str = f"{float(f_corr):.{max_decimals}f}"

        return f_ref_str, f_corr_str


    def show_parameters(self):
        """Collect and display all parameters when Validate is clicked."""

        # Ensure that CBS F values are valid if CBS is selected
        if self.cbs_checkbox.isChecked() and self.cbs_option.currentText() == "Other":
            if not self.validate_f_values_immediate():
                return  # Stop if F values are not valid

        # Check if job name is valid
        if not self.job_name_input.text().strip():
            QMessageBox.warning(self, "Missing Job Name", "Please assign a valid job name before proceeding.")
            return

        # Set the job name
        self.job_name = self.job_name_input.text().strip()

        # Ensure that `ledaw_out_root` is not overwritten if the user manually set it
        if not self.is_manual_directory:  # Only set if the user didn't manually update the directory
            self.ledaw_out_root = rf"./LEDAW-OUT/{self.job_name}"

        # Update the output directory text box with the LEDAW output root
        self.output_dir_input.blockSignals(True)  # Prevent feedback loop when setting text
        self.output_dir_input.setText(self.ledaw_out_root)
        self.output_dir_input.blockSignals(False)

        # Ensure that at least one of N-body, Two-body, and plot LED is selected
        if not self.check_selections():
            return 

        # Collect method, CBS, CPS, and other parameter values
        method = self.validate_method()
        if method is None:
            return  # Do not proceed if method is invalid

        conversion_factor = self.conversion_factor
        self.use_ref_as_rhf_in_hfld = self.get_rhf_selection()

        F = self.F if hasattr(self, 'F') else None  # Get CPS F value
        
        F_ref_cbs, F_corr_cbs = self.validate_f_values()
        # Format F values for display
        if F_ref_cbs is not None and F_corr_cbs is not None:
            f_ref_display, f_corr_display = self.format_f_values_for_display(str(F_ref_cbs), str(F_corr_cbs))
        else:
            f_ref_display, f_corr_display = None, None

        # Checkboxes to determine what to perform
        perform_nbody = self.nbody_checkbox.isChecked()
        perform_twobody = self.twobody_checkbox.isChecked()
        perform_cps = self.cps_checkbox.isChecked()
        perform_cbs = self.cbs_checkbox.isChecked()
        perform_plot = self.plot_checkbox.isChecked()

        # Create the confirmation dialog window
        confirm_dialog = QDialog(self)
        confirm_dialog.setWindowTitle("Confirm Parameters")
        confirm_dialog.setMinimumSize(500, 400)

        # Remove the "?" help button
        confirm_dialog.setWindowFlags(confirm_dialog.windowFlags() & ~Qt.WindowType.WindowContextHelpButtonHint)

        # Create a layout for the dialog
        layout = QVBoxLayout()

        # Add QLabel with HTML formatting for bold text
        label = QLabel("If you want to continue with the following parameters?<br>"
                       "Please press OK. The Home tab will then be locked to proceed.")
        layout.addWidget(label)

        # Warn for nonempty directories
        warning_message = ""
        if os.path.exists(self.ledaw_out_root) and os.listdir(self.ledaw_out_root):
            warning_message = "<br><b>Warning:</b> The folder specified for writing LEDAW outputs is not empty. To prevent interference and avoid overwriting existing files, it is strongly recommended to press 'Cancel' and provide an empty directory path.<br><br>"
        
        # Create a QTextEdit for displaying the parameters, set it to read-only and allow selectable text
        text_edit = QTextEdit()
        text_edit.setReadOnly(True)
        text_edit.setHtml(f"""{warning_message}
            <b>LEDAW Output Directory:</b> {self.ledaw_out_root}<br>
            <b>Method:</b> {method}<br>
            <b>Conversion Factor:</b> {conversion_factor}<br>
            <b>Use RHF Reference for Subsystems in HFLD:</b> {self.use_ref_as_rhf_in_hfld}<br>
            <b>Perform N-body LED:</b> {perform_nbody}<br>
            <b>Perform Two-body LED:</b> {perform_twobody}<br>
            <b>Perform CPS Extrapolation:</b> {perform_cps}<br>
            <b>F for CPS:</b> {F}<br>
            <b>Perform CBS Extrapolation:</b> {perform_cbs}<br>
            <b>F for CBS (Reference Part):</b> {f_ref_display}<br>
            <b>F for CBS (Correlation Part):</b> {f_corr_display}<br>
            <b>Perform Plot LED Heat Maps:</b> {perform_plot}<br>""")

        layout.addWidget(text_edit)

        # Create OK and Cancel buttons
        button_layout = QHBoxLayout()
        ok_button = QPushButton("OK")
        cancel_button = QPushButton("Cancel")
        button_layout.addWidget(ok_button)
        button_layout.addWidget(cancel_button)
        layout.addLayout(button_layout)

        # Connect buttons
        ok_button.clicked.connect(lambda: self.handle_confirmation(confirm_dialog))
        cancel_button.clicked.connect(confirm_dialog.reject)

        # Set the layout for the dialog
        confirm_dialog.setLayout(layout)

        # Show the dialog
        if confirm_dialog.exec() == QDialog.DialogCode.Accepted:
            # Add tabs in order based on user selection
            self.add_tabs_in_order()
        else:
            # If No, return to the Home tab
            self.tabs.setCurrentIndex(0)


    def handle_confirmation(self, dialog):
        """Lock Home tab and delete tmp directory created in previous run at the same directory"""
        dialog.accept()  # Close the dialog


    def show_nbody_info(self):
        msg = QMessageBox()
        msg.setWindowTitle("N-body LED Info")

        try:
            # If running under PyInstaller, use sys._MEIPASS for the path
            if getattr(sys, 'frozen', False):  # Running as a PyInstaller bundle
                base_path = sys._MEIPASS
            else:
                base_path = os.path.dirname(os.path.abspath(sys.argv[0]))

            # Set the icon based on the operating system
            if sys.platform.startswith('win'):
                icon_path = os.path.join(base_path, 'ledaw_package', 'favicon.ico')
            elif sys.platform == 'darwin':  # macOS
                icon_path = os.path.join(base_path, 'ledaw_package', 'favicon.icns')
            else:  # Assume Linux or other OS
                icon_path = os.path.join(base_path, 'ledaw_package', 'favicon.png')

            # Set the message box icon
            msg.setWindowIcon(QIcon(icon_path))

        except Exception as e:
            print(f"Could not set the message box icon: {e}")
            # Continue without crashing if the icon fails to load
    
        msg.setText("To perform this analysis, ORCA output files are required for both the supersystem (i.e., the adduct consisting of two or more fragments, on which interaction energy is calculated) and the subsystems used to compute the interaction energy. No LED data are needed in the ORCA output files for subsystems that consist of only one fragment.")
        msg.exec()


    def show_twobody_info(self):
        msg = QMessageBox()
        msg.setWindowTitle("Two-body LED Info")

        try:
            # If running under PyInstaller, use sys._MEIPASS for the path
            if getattr(sys, 'frozen', False):  # Running as a PyInstaller bundle
                base_path = sys._MEIPASS
            else:
                base_path = os.path.dirname(os.path.abspath(sys.argv[0]))

            # Set the icon based on the operating system
            if sys.platform.startswith('win'):
                icon_path = os.path.join(base_path, 'ledaw_package', 'favicon.ico')
            elif sys.platform == 'darwin':  # macOS
                icon_path = os.path.join(base_path, 'ledaw_package', 'favicon.icns')
            else:  # Assume Linux or other OS
                icon_path = os.path.join(base_path, 'ledaw_package', 'favicon.png')

            # Set the message box icon
            msg.setWindowIcon(QIcon(icon_path))

        except Exception as e:
            print(f"Could not set the message box icon: {e}")
            # Continue without crashing if the icon fails to load
    
        msg.setText("To perform this analysis, ORCA output files are required for all the two-fragment subsystems of a supersystem, whose fragments belong to different subsystems, and for all the corresponding monomers. No LED data are needed in the ORCA output files of monomers. The two-body module collects LED interaction energy components of the isolated two-fragment subsystems on NxN matrices.\n\n"
                   "If both N-body and two-body analyses are selected, cooperativity analysis will be performed as well.")
        msg.exec()


    def show_cps_info(self):
        msg = QMessageBox()
        msg.setWindowTitle("CPS Extrapolation Info")

        try:
            # If running under PyInstaller, use sys._MEIPASS for the path
            if getattr(sys, 'frozen', False):  # Running as a PyInstaller bundle
                base_path = sys._MEIPASS
            else:
                base_path = os.path.dirname(os.path.abspath(sys.argv[0]))

            # Set the icon based on the operating system
            if sys.platform.startswith('win'):
                icon_path = os.path.join(base_path, 'ledaw_package', 'favicon.ico')
            elif sys.platform == 'darwin':  # macOS
                icon_path = os.path.join(base_path, 'ledaw_package', 'favicon.icns')
            else:  # Assume Linux or other OS
                icon_path = os.path.join(base_path, 'ledaw_package', 'favicon.png')

            # Set the message box icon
            msg.setWindowIcon(QIcon(icon_path))

        except Exception as e:
            print(f"Could not set the message box icon: {e}")
            # Continue without crashing if the icon fails to load
    
        msg.setText("To perform CPS extrapolation, ORCA output files containing results with both a looser PNO and "
                    "a tighter PNO settings are necessary, with adjacent TCutPNO exponents, such as, TCutPNO = 1e-6 and 1e-7, respectively.")
        msg.exec()


    def show_cbs_info(self):
        msg = QMessageBox()
        msg.setWindowTitle("CBS Extrapolation Info")

        try:
            # If running under PyInstaller, use sys._MEIPASS for the path
            if getattr(sys, 'frozen', False):  # Running as a PyInstaller bundle
                base_path = sys._MEIPASS
            else:
                base_path = os.path.dirname(os.path.abspath(sys.argv[0]))

            # Set the icon based on the operating system
            if sys.platform.startswith('win'):
                icon_path = os.path.join(base_path, 'ledaw_package', 'favicon.ico')
            elif sys.platform == 'darwin':  # macOS
                icon_path = os.path.join(base_path, 'ledaw_package', 'favicon.icns')
            else:  # Assume Linux or other OS
                icon_path = os.path.join(base_path, 'ledaw_package', 'favicon.png')

            # Set the message box icon
            msg.setWindowIcon(QIcon(icon_path))

        except Exception as e:
            print(f"Could not set the message box icon: {e}")
            # Continue without crashing if the icon fails to load
    
        msg.setText("To perform CBS extrapolation, ORCA output files containing results with both a smaller basis set "
                    "and a larger basis set with adjacent cardinal numbers are necessary, such as, cc-pVTZ and cc-pVQZ, respectively.\n\n"
                    "For CBS(2/3) and (3/4), F coefficients are chosen as in the Supporting Information of the open-shell HFLD method paper "
                    "(DOI: 10.1021/acs.jctc.1c01295).")
        msg.exec()


    def show_plot_info(self):
        msg = QMessageBox()
        msg.setWindowTitle("Plot LED Heat Map Info")

        try:
            # If running under PyInstaller, use sys._MEIPASS for the path
            if getattr(sys, 'frozen', False):  # Running as a PyInstaller bundle
                base_path = sys._MEIPASS
            else:
                base_path = os.path.dirname(os.path.abspath(sys.argv[0]))

            # Set the icon based on the operating system
            if sys.platform.startswith('win'):
                icon_path = os.path.join(base_path, 'ledaw_package', 'favicon.ico')
            elif sys.platform == 'darwin':  # macOS
                icon_path = os.path.join(base_path, 'ledaw_package', 'favicon.icns')
            else:  # Assume Linux or other OS
                icon_path = os.path.join(base_path, 'ledaw_package', 'favicon.png')

            # Set the message box icon
            msg.setWindowIcon(QIcon(icon_path))

        except Exception as e:
            print(f"Could not set the message box icon: {e}")
            # Continue without crashing if the icon fails to load
    
        msg.setText("It is necessary to select one of N-body LED and two-body LED or to have files of these analyses from previous runs in the LEDAW output directory of the given job name.")
        msg.exec()


    def create_manage_session_tab(self):
        """Create the Manage Session tab with Resume, Start Over, and Exit options."""
        manage_session_tab = QWidget()
        layout = QVBoxLayout()

        # Create Continue, Start Over, and Exit buttons
        continue_button = QPushButton("Resume")
        continue_button.setFixedSize(100, 40)
        continue_button.clicked.connect(self.resume_action)
        continue_button.setStyleSheet("QPushButton {background-color: #ebedef;} QPushButton:hover {background-color: paleturquoise;}")

        restart_button = QPushButton("Start Over")
        restart_button.setFixedSize(100, 40)
        restart_button.clicked.connect(self.reset_application)
        restart_button.setStyleSheet("QPushButton {background-color: #ebedef;} QPushButton:hover {background-color: paleturquoise;}")
        
        exit_button = QPushButton("Exit")
        exit_button.setFixedSize(100, 40)
        exit_button.clicked.connect(self.confirm_exit)
        exit_button.setStyleSheet("QPushButton {background-color: #ebedef;} QPushButton:hover {background-color: paleturquoise;}")

        # Center buttons in the middle of the tab
        button_layout = QHBoxLayout()
        button_layout.addStretch()  # Add space before buttons to center them
        button_layout.addWidget(continue_button)
        button_layout.addWidget(restart_button)
        button_layout.addWidget(exit_button)
        button_layout.addStretch()  # Add space after buttons to center them

        # Add buttons to the Manage Session tab layout
        layout.addStretch()  # Add stretch at the top to push the buttons to the center vertically
        layout.addLayout(button_layout)
        layout.addStretch()  # Add stretch at the bottom to keep the buttons centered

        manage_session_tab.setLayout(layout)
        self.tabs.addTab(manage_session_tab, "Manage Session")


    def resume_action(self):
        """Action for Resume button: Go to the last dynamically created tab."""
        if self.created_tabs:
            # Find the lowest index of the dynamically created tabs
            highest_index = max(self.tabs.indexOf(tab) for tab in self.created_tabs)

            # Set the current tab to the first dynamically created tab
            self.tabs.setCurrentIndex(highest_index)
        else:
            # If no dynamic tabs exist, fallback to Home tab
            self.tabs.setCurrentIndex(self.tabs.indexOf(self.home_tab))


    def reset_application(self):
        """Reset the application to its initial state."""
        self.bypass_validation = True  # Enable bypass validation during reset

        # Remove all dynamically created tabs (N-body, Two-body, Plot)
        self.remove_existing_tabs()

        self.nbody_locked = False
        self.twobody_locked = False 
        self.plot_locked = False 

        # Clear the list of created tabs
        self.created_tabs.clear()

        # Reset all selections and variables
        self.job_name_input.clear()
        self.output_dir_input.clear()
        self.method.setCurrentIndex(0)
        self.unit_choice.setCurrentIndex(0)
        self.nbody_checkbox.setChecked(False)
        self.twobody_checkbox.setChecked(False)
        self.cps_checkbox.setChecked(False)
        self.cbs_checkbox.setChecked(False)
        self.plot_checkbox.setChecked(False)
        self.use_rhf_checkbox.setChecked(False)
        self.last_selected_dir = './'
        self.plot_unlock_count = 0
        
        if hasattr(self, 'show_diag_checkbox'):
            self.show_diag_checkbox.setChecked(False)

        if hasattr(self, 'delete_old_plots_checkbox'):
            self.delete_old_plots_checkbox.setChecked(False)
            
        if hasattr(self, 'rerun_plot_copy_button'):
            # Remove the copy button if it exists
            parent_layout = self.rerun_plot_button.parentWidget().layout()
            parent_layout.removeWidget(self.rerun_plot_copy_button)
            self.rerun_plot_copy_button.deleteLater()  # Clean up the copy button
            del self.rerun_plot_copy_button  # Remove reference to the copy button
        
        # Reset CBS-specific fields
        self.cbs_option.setCurrentIndex(0)  # Reset CBS option to default (CBS(3/4))
        self.f_ref_input.clear()  # Clear F ref input
        self.f_corr_input.clear()  # Clear F corr input
        self.set_cbs_defaults()  # Reapply CBS defaults based on the reset selection
        self.f_ref = None
        self.f_corr = None
        self.F = None
        
        # Hide and reset CBS and RHF options
        self.cbs_label.setVisible(False)
        self.cbs_option.setVisible(False)
        self.f_ref_input.setVisible(False)
        self.f_corr_input.setVisible(False)

        # Reset relabel_mapping explicitly
        self.relabel_mappings = {}

        # Re-enable the Home tab elements
        self.reset_home_tab()  # <-- Reset the Home tab to be editable again
        self.tabs.setCurrentIndex(self.tabs.indexOf(self.home_tab))  # Return to the Home tab
        
        # Delete the tmp directory if it exists
        self.delete_tmp_files()

        self.bypass_validation = False  # Disable bypass after reset


    def remove_existing_tabs(self):
        """Remove any existing N-body, Two-body, and Plot tabs to avoid duplication."""
        for tab in self.created_tabs:
            index = self.tabs.indexOf(tab)
            if index != -1:
                self.tabs.removeTab(index)

        # Clear the list of created tabs once all are removed
        self.created_tabs.clear()

        # Also reset the tab references to None so they can be recreated
        self.nbody_tab = None
        self.twobody_tab = None
        self.plot_tab = None
        self.run_tab = None


    def create_nbody_tab(self):
        """Create the N-body tab with different content based on CPS and CBS selections."""
        # Check CPS and CBS selections
        cps_selected = self.cps_checkbox.isChecked()
        cbs_selected = self.cbs_checkbox.isChecked()

        # Create the N-body tab widget
        self.nbody_tab = QWidget()

        # Create a scroll area for the layout
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)  # Allow the scroll area to resize with the window

        # Create a widget that will be placed inside the scroll area
        scroll_widget = QWidget()

        # Create a layout for the scroll widget
        layout = QVBoxLayout()
        scroll_widget.setLayout(layout)

        # Scenario 1: Both CPS and CBS selected
        if cps_selected and cbs_selected:
            self.add_cps_cbs_nbody_layout(layout)

        # Scenario 2: Only CPS selected
        elif cps_selected and not cbs_selected:
            self.add_cps_nbody_layout(layout)

        # Scenario 3: Only CBS selected
        elif not cps_selected and cbs_selected:
            self.add_cbs_nbody_layout(layout)

        # Scenario 4: Neither CPS nor CBS selected
        else:
            # Call the function that contains the logic for the standard N-body tab (without CPS or CBS)
            self.add_standard_nbody_layout(layout)

        # Set the scroll widget as the content of the scroll area
        scroll_area.setWidget(scroll_widget)

        # Now, add the scroll area to the N-body tab instead of setting the layout directly
        layout_wrapper = QVBoxLayout(self.nbody_tab)
        layout_wrapper.addWidget(scroll_area)

        # Insert the N-body tab before Exit and About
        exit_index = self.tabs.count() - 2  # Assuming Exit is the second last tab
        self.tabs.insertTab(exit_index, self.nbody_tab, "N-body")
        self.created_tabs.append(self.nbody_tab)

        # Move to the N-body tab after Home confirmation
        self.tabs.setCurrentIndex(self.tabs.indexOf(self.nbody_tab))

        # Return the correct tab index
        return self.tabs.indexOf(self.nbody_tab)


    def ensure_list_size(self, file_list, index):
        """Ensure the list is large enough to accommodate a new file at the given index."""
        while len(file_list) <= index:
            file_list.append('')


    def filter_empty_rows_and_validate_main(self, main_files, alt_files, section):
        """Filter out empty rows from filenames and check if alt file is provided without a main file.
           Ensure the supersystem field (index 0) is preserved even if cleared, but prompt user on confirm if it's empty.
        """
        filtered_main = []
        filtered_alt = []

        # If the supersystem main file is empty, preserve it as '' but prompt the user during confirmation
        if not main_files[0].strip():  # If the supersystem main file is empty
            QMessageBox.warning(self, "Missing Supersystem File",
                                f"Please specify files for all supersystems.")
            return None, None  # Stop processing and return early if invalid

        # Ensure that the alternative supersystem file (index 0) can be empty, but not alt without main
        if alt_files[0] and not main_files[0]:
            QMessageBox.warning(self, "Missing Main File",
                                f"You have provided an alternative file without a corresponding supersystem main file for section {section}.")
            return None, None  # Stop processing and return early if invalid

        # Add the supersystem (index 0) without modification (even if it's an empty string '')
        filtered_main.append(main_files[0])
        filtered_alt.append(alt_files[0])

        # Loop through the rest of the subsystem files, starting from index 1
        for main, alt in zip(main_files[1:], alt_files[1:]):
            if alt and not main:  # Alt file is present but no main file
                QMessageBox.warning(self, "Missing Main File",
                                    "You have provided an alternative file without a corresponding main file. Please ensure that each subsystem has a main file.")
                return None, None  # Stop processing if invalid
            if main or alt:  # If at least one of the fields is filled, it's a valid entry
                filtered_main.append(main)
                filtered_alt.append(alt)

        return filtered_main, filtered_alt


    def count_valid_subsystems(self, main_files):
        """Count the number of non-empty subsystem files, skipping the first supersystem entry."""
        return len([f for f in main_files[1:] if f])  # Skip the first file (supersystem) and count non-empty subsystem files


    def show_relabel_info(self):
        msg = QMessageBox()
        msg.setWindowTitle("Fragment Reordering Info")

        try:
            # If running under PyInstaller, use sys._MEIPASS for the path
            if getattr(sys, 'frozen', False):  # Running as a PyInstaller bundle
                base_path = sys._MEIPASS
            else:
                base_path = os.path.dirname(os.path.abspath(sys.argv[0]))

            # Set the icon based on the operating system
            if sys.platform.startswith('win'):
                icon_path = os.path.join(base_path, 'ledaw_package', 'favicon.ico')
            elif sys.platform == 'darwin':  # macOS
                icon_path = os.path.join(base_path, 'ledaw_package', 'favicon.icns')
            else:  # Assume Linux or other OS
                icon_path = os.path.join(base_path, 'ledaw_package', 'favicon.png')

            # Set the message box icon
            msg.setWindowIcon(QIcon(icon_path))

        except Exception as e:
            print(f"Could not set the message box icon: {e}")
            # Continue without crashing if the icon fails to load
    
        msg.setText("Note 1: Relabel mapping list is used to reorder fragment labels in the main supersystem (adduct) ORCA output file if the user is not satisfied with the original fragment labels in this file.\n\n"
                    "Example: 5,3,2,1,4 \nIn this example, fragment originally labeled as 1 in the main supersystem ORCA file will be relabeled as 5; 2 as 3; 3 as 2; 4 as 1; and 5 as 4.\n\n"
                    "Note 2: The list must contain a complete sequence of integers from 1 up to the highest fragment label present in the main supersystem file, without any missing integers. Note that if fragments are not "
					"labeled consecutively from 1 to the total number of fragments in the supersystem file, the highest label will differ from the actual fragment count. The LEDAW GUI only checks whether the list includes "
					"consecutive integers starting from 1, but it does not verify the highest label. Therefore, if the length of relabel mapping list is mistakenly provided differently than the highest fragment label, this "
					"will interrupt LEDAW’s execution with an error.\n\n"					
					"Note 3: Regardless of whether a relabel mapping is used, users do not need to worry if the same fragment has different labels across the supersystem main, supersystem alternative, and subsystem ORCA output files. "
                    "All fragment labels are automatically standardized to match the main supersystem file - or its reordered labels if a relabel mapping is applied.")
        msg.exec()


    def add_relabel_mapping_section(self, target_layout, section):
        """Adds a relabel mapping section that dynamically shows/hides the input field when the checkbox is toggled."""

        if not hasattr(self, 'relabel_mappings'):
            self.relabel_mappings = {}  # Dictionary to store relabel mappings for each section
        if not hasattr(self, 'relabel_inputs'):
            self.relabel_inputs = {}  # Dictionary to store relabel input fields for each section
        self.relabel_mappings[section] = None  # Initially set to None for the section

        # Fragment Reordering Checkbox
        relabel_checkbox = QCheckBox("Do you want to reorder fragment labels in the main supersystem file?")

        # Info Button for the relabel section
        relabel_info_button = QPushButton('i')
        relabel_info_button.setFixedSize(25, 25)  # Minimal size
        relabel_info_button.setStyleSheet("border-radius: 12px; background-color: lightgray; color: black; font-weight: bold; font-size: 13px;")

        # Connect the info button to show the relabel information
        relabel_info_button.clicked.connect(self.show_relabel_info)

        # Layout for the checkbox and info button on the same row
        relabel_checkbox_layout = QHBoxLayout()
        relabel_checkbox_layout.addWidget(relabel_checkbox)
        relabel_checkbox_layout.addStretch()  # Push the info button to the right
        relabel_checkbox_layout.addWidget(relabel_info_button)

        # Add the checkbox and info button row to the layout
        target_layout.addLayout(relabel_checkbox_layout)

        # Create a horizontal layout for the relabel input field (it will occupy the entire row)
        relabel_input_row_layout = QHBoxLayout()

        # Create the relabel input field
        relabel_input = QLineEdit()
        relabel_input.setPlaceholderText("Enter relabel mapping list")

        # Set the size policy to ensure it occupies the entire row and resizes dynamically
        relabel_input.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        relabel_input.setFixedHeight(30)  # Adjust height to match other UI components

        # Add the input field to the row layout
        relabel_input_row_layout.addWidget(relabel_input)

        # Add the row layout to the main layout (target_layout)
        target_layout.addLayout(relabel_input_row_layout)

        # Store the input field for this section
        self.relabel_inputs[section] = relabel_input

        # Initially hide the input field, but keep the layout intact
        relabel_input.setVisible(False)

        
        # Function to validate relabel mapping input after typing is finished
        def validate_relabel_mapping():
            relabel_input = self.relabel_inputs.get(section)  # Get the correct input field for the section
            if relabel_input is None:
                return None

            text = relabel_input.text().strip()  # Use the section-specific input
            if not text:  # If the input is empty, treat it as None
                self.relabel_mappings[section] = None
                return None

            try:
                # Split and parse the input as integers
                relabel_mapping = [int(x) for x in text.strip('[]').split(',')]

                # Check for duplicates
                if len(relabel_mapping) != len(set(relabel_mapping)):
                    raise ValueError("Duplicate values found in the relabel mapping.")

                # Check that numbers range from 1 to n with no missing values
                if min(relabel_mapping) != 1 or max(relabel_mapping) != len(relabel_mapping):
                    raise ValueError("Relabel mapping must be a complete sequence from 1 to the highest fragment label—no numbers can be missing.")

                # If all checks pass, store the relabel mapping
                self.relabel_mappings[section] = relabel_mapping
                return relabel_mapping

            except ValueError as e:
                # Show the specific error message based on the raised exception
                QMessageBox.warning(self, "Invalid Input", f"{str(e)} Please correct it. Otherwise, relabel mapping will not be applied.")
                self.relabel_mappings[section] = None  # Reset mapping if invalid
                return None

            
        # Function to handle the checkbox toggle event (show/hide the input field)
        def toggle_relabel_input(checked):
            relabel_input.setVisible(checked)  # Show the QLineEdit if checked, hide it if unchecked
            if not checked:
                # Clear the input when hidden and reset the validation
                self.relabel_mappings[section] = None
            else:
                # If the checkbox is checked, trigger validation when the user finishes editing the field
                relabel_input.editingFinished.connect(validate_relabel_mapping)

        # Connect the checkbox toggle event to show/hide the input field
        relabel_checkbox.toggled.connect(toggle_relabel_input)


    def add_cps_cbs_nbody_layout(self, layout):
        """Create a CPS and CBS layout for the N-body tab"""

        self.nbody_tab.setMinimumSize(700, 400)

        # Initialize lists to hold file paths for CPS and CBS sections
        self.main_filenames_SB_LPNO = ['']
        self.alternative_filenames_SB_LPNO = ['']
        self.main_filenames_SB_TPNO = ['']
        self.alternative_filenames_SB_TPNO = ['']
        self.main_filenames_LB_LPNO = [""] 
        self.alternative_filenames_LB_LPNO = ['']
        self.main_filenames_LB_TPNO = [""] 
        self.alternative_filenames_LB_TPNO = ['']

        # Info button to be placed top-right
        info_button = QPushButton('i')
        info_button.setFixedSize(25, 25)
        info_button.setStyleSheet("border-radius: 12px; background-color: lightgray; color: black; font-weight: bold; font-size: 13px;")


        def show_info():
            QMessageBox.information(self, "Alternative ORCA Output File Info",
                                    "If part of the LED data is missing in an ORCA output file but present in another ORCA output file for the same system, both files can be specified as complementary to each other.\n\n"
                                    "For example, if DLPNO-CCSD(T) is chosen, the DLPNO-CCSD output can be used as the main file, while an alternative file that contains the triples contribution but lacks the HF/LED decomposition is specified. Similarly, an HFLD/LED output can be provided as the alternative file for the HF/LED portion.\n\n"
                                    "Note 1: If data for an LED component is present in both the main and alternative files, the main file will be prioritized.\n\n"
                                    "Note 2: As an example, if you perform the calculations with cc-pvTZ and cc-pVQZ basis sets and with TCutPNO= 1e-6 and TCutPNO=1e-7, then:\n\n"
                                    "Smaller Basis Set: cc-pVTZ ; Larger Basis Set: cc-pVQZ\n"
                                    "Looser PNO: TCutPNO = 1e-6 ; Tighter PNO: TCutPNO = 1e-7")

        info_button.clicked.connect(show_info)

        # Add the info button layout above the header layout
        info_layout = QHBoxLayout()
        info_layout.addStretch()  # Push the info button to the right
        info_layout.addWidget(info_button)
        layout.addLayout(info_layout)  # Info button on its own row

        # Spacer between Main and Alternative headers with a fixed size of 3px
        header_spacer = QSpacerItem(3, 0, QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Minimum)


        def select_file(button, file_type, row, section, is_supersystem=False):
            file_path, _ = QFileDialog.getOpenFileName(self, "Select ORCA Output File", self.last_selected_dir)

            # Check if the user pressed cancel (file_path will be empty)
            if not file_path:
                return  # Do nothing if no file was selected (i.e., cancel was pressed)

            # Normalize the selected file path
            file_path = self.normalize_path_gui(file_path)

            # Update the last selected directory
            self.last_selected_dir = file_path.rsplit("/", 1)[0]

            # Use the ensure_list_size method from the class
            if section == "SB_LPNO":
                self.ensure_list_size(self.main_filenames_SB_LPNO, row)
                self.ensure_list_size(self.alternative_filenames_SB_LPNO, row)
                if file_type == "main":
                    self.main_filenames_SB_LPNO[row] = file_path
                else:
                    self.alternative_filenames_SB_LPNO[row] = file_path or ''

            elif section == "SB_TPNO":
                self.ensure_list_size(self.main_filenames_SB_TPNO, row)
                self.ensure_list_size(self.alternative_filenames_SB_TPNO, row)
                if file_type == "main":
                    self.main_filenames_SB_TPNO[row] = file_path
                else:
                    self.alternative_filenames_SB_TPNO[row] = file_path or ''

            elif section == "LB_LPNO":
                self.ensure_list_size(self.main_filenames_LB_LPNO, row)
                self.ensure_list_size(self.alternative_filenames_LB_LPNO, row)
                if file_type == "main":
                    self.main_filenames_LB_LPNO[row] = file_path
                else:
                    self.alternative_filenames_LB_LPNO[row] = file_path or ''

            elif section == "LB_TPNO":
                self.ensure_list_size(self.main_filenames_LB_TPNO, row)
                self.ensure_list_size(self.alternative_filenames_LB_TPNO, row)
                if file_type == "main":
                    self.main_filenames_LB_TPNO[row] = file_path
                else:
                    self.alternative_filenames_LB_TPNO[row] = file_path or ''

            # Update the button text to show the selected file path
            button.setText(file_path)


        def add_subsystems(target_layout, section):
            """Add and process subsystem file paths"""

            def select_multiple_files(section, current_index):
                """Select multiple files and create rows dynamically."""
                file_paths, _ = QFileDialog.getOpenFileNames(self, "Select ORCA Output Files", self.last_selected_dir)

                # Check if the user pressed cancel (file_paths will be empty)
                if not file_paths:
                    return  # Do nothing if no file was selected (i.e., cancel was pressed)

                # Normalize all selected file paths
                file_paths = [self.normalize_path_gui(path) for path in file_paths]

                # Update the last selected directory
                self.last_selected_dir = file_paths[0].rsplit("/", 1)[0] if file_paths else self.last_selected_dir

                # Process each selected file as if the "Add Subsystem Files" button was pressed for each
                for file_path in file_paths:
                    add_single_subsystem_row(file_path, section, current_index)
                    current_index += 1


            def add_single_subsystem_row(file_path, section, row):
                """Create a single row for a subsystem with main and alternative file options."""
                subsystem_row = QHBoxLayout()

                subsystem_main = QLineEdit()
                subsystem_alt = QLineEdit()

                subsystem_main.setPlaceholderText("Subsystem Main File")
                subsystem_alt.setPlaceholderText("Subsystem Alternative File")

                subsystem_main.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
                subsystem_alt.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
                subsystem_main.setFixedHeight(30)
                subsystem_alt.setFixedHeight(30)

                subsystem_main_button = QPushButton("Select Subsystem File")
                subsystem_main_button.setFixedWidth(150)
                subsystem_main_button.setFixedHeight(30)

                subsystem_alt_button = QPushButton("Select Subsystem File")
                subsystem_alt_button.setFixedWidth(150)
                subsystem_alt_button.setFixedHeight(30)

                # Automatically populate the main file with the selected file path
                subsystem_main.setText(file_path)

                if section == "SB_LPNO":
                    self.ensure_list_size(self.main_filenames_SB_LPNO, row)
                    self.ensure_list_size(self.alternative_filenames_SB_LPNO, row)
                    self.main_filenames_SB_LPNO[row] = file_path

                    subsystem_main_button.clicked.connect(lambda: select_file(subsystem_main, "main", row, section))
                    subsystem_alt_button.clicked.connect(lambda: select_file(subsystem_alt, "alt", row, section))

                    subsystem_main.editingFinished.connect(lambda: self.update_filename(self.main_filenames_SB_LPNO, subsystem_main.text(), row))
                    subsystem_alt.editingFinished.connect(lambda: self.update_filename(self.alternative_filenames_SB_LPNO, subsystem_alt.text(), row))

                elif section == "SB_TPNO":
                    self.ensure_list_size(self.main_filenames_SB_TPNO, row)
                    self.ensure_list_size(self.alternative_filenames_SB_TPNO, row)
                    self.main_filenames_SB_TPNO[row] = file_path

                    subsystem_main_button.clicked.connect(lambda: select_file(subsystem_main, "main", row, section))
                    subsystem_alt_button.clicked.connect(lambda: select_file(subsystem_alt, "alt", row, section))

                    subsystem_main.editingFinished.connect(lambda: self.update_filename(self.main_filenames_SB_TPNO, subsystem_main.text(), row))
                    subsystem_alt.editingFinished.connect(lambda: self.update_filename(self.alternative_filenames_SB_TPNO, subsystem_alt.text(), row))

                elif section == "LB_LPNO":
                    self.ensure_list_size(self.main_filenames_LB_LPNO, row)
                    self.ensure_list_size(self.alternative_filenames_LB_LPNO, row)
                    self.main_filenames_LB_LPNO[row] = file_path

                    subsystem_main_button.clicked.connect(lambda: select_file(subsystem_main, "main", row, section))
                    subsystem_alt_button.clicked.connect(lambda: select_file(subsystem_alt, "alt", row, section))

                    subsystem_main.editingFinished.connect(lambda: self.update_filename(self.main_filenames_LB_LPNO, subsystem_main.text(), row))
                    subsystem_alt.editingFinished.connect(lambda: self.update_filename(self.alternative_filenames_LB_LPNO, subsystem_alt.text(), row))

                elif section == "LB_TPNO":
                    self.ensure_list_size(self.main_filenames_LB_TPNO, row)
                    self.ensure_list_size(self.alternative_filenames_LB_TPNO, row)
                    self.main_filenames_LB_TPNO[row] = file_path

                    subsystem_main_button.clicked.connect(lambda: select_file(subsystem_main, "main", row, section))
                    subsystem_alt_button.clicked.connect(lambda: select_file(subsystem_alt, "alt", row, section))

                    subsystem_main.editingFinished.connect(lambda: self.update_filename(self.main_filenames_LB_TPNO, subsystem_main.text(), row))
                    subsystem_alt.editingFinished.connect(lambda: self.update_filename(self.alternative_filenames_LB_TPNO, subsystem_alt.text(), row))

                # Add the buttons and text fields to the layout
                subsystem_row.addWidget(subsystem_main_button)
                subsystem_row.addWidget(subsystem_main)
                subsystem_row.addWidget(subsystem_alt_button)
                subsystem_row.addWidget(subsystem_alt)

                # Add the row to the target layout
                target_layout.addLayout(subsystem_row)

            # Get the current number of subsystems to maintain indexing consistency
            current_index = len(self.main_filenames_SB_LPNO) if section == "SB_LPNO" else \
                            len(self.main_filenames_SB_TPNO) if section == "SB_TPNO" else \
                            len(self.main_filenames_LB_LPNO) if section == "LB_LPNO" else \
                            len(self.main_filenames_LB_TPNO)

            # Trigger the file selection and row creation process
            select_multiple_files(section, current_index)

        # Adding the "Smaller Basis Set and Looser PNO" Section
        sb_lpno_label = QLabel("Smaller Basis Set and Looser PNO")
        sb_lpno_label.setStyleSheet("background-color: #d3e2f4; font-size: 12px; padding: 5px; text-align: center; font-weight: bold;")
        sb_lpno_label.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        layout.addWidget(sb_lpno_label)

        # Add relabel mapping section right after the heading
        self.add_relabel_mapping_section(layout, "SB_LPNO")

        # Create header labels for main and alternative sections
        main_label = QLabel("Main ORCA Output Files")
        alt_label = QLabel("Alternative ORCA Output Files (Optional)")

        # Styling the labels: equal size, bold text (font 2px smaller than select button text), and height matching select boxes
        main_label.setStyleSheet("background-color: #ebedef; font-size: 11px; padding: 5px; text-align: center; font-weight: bold;")
        alt_label.setStyleSheet("background-color: #ebedef; font-size: 11px; padding: 5px; text-align: center; font-weight: bold;")

        # Set equal size policy for dynamic resizing of both headers
        main_label.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        alt_label.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)

        # Match the height of the select boxes and headings (assuming 30px for the select buttons)
        select_button_height = 30
        main_label.setFixedHeight(select_button_height)
        alt_label.setFixedHeight(select_button_height)

        # Layout for the headers
        sb_lpno_header_layout = QHBoxLayout()
        sb_lpno_header_layout.addWidget(main_label)  # Main label dynamically resizes
        sb_lpno_header_layout.addSpacerItem(header_spacer)  # Fixed 3px spacer between labels
        sb_lpno_header_layout.addWidget(alt_label)  # Alt label dynamically resizes

        # Add the header layout to the main layout (below the info button)
        layout.addLayout(sb_lpno_header_layout)

        # Supersystem Section (Fixed-size select buttons, expandable populate fields)
        sb_lpno_supersystem_layout = QHBoxLayout()

        self.sb_lpno_supersystem_main = QLineEdit()
        self.sb_lpno_supersystem_alt = QLineEdit()

        # Set placeholders
        self.sb_lpno_supersystem_main.setPlaceholderText("Supersystem Main File")
        self.sb_lpno_supersystem_alt.setPlaceholderText("Supersystem Alternative File")

        # Set size policy and height to match the heading size
        self.sb_lpno_supersystem_main.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.sb_lpno_supersystem_alt.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.sb_lpno_supersystem_main.setFixedHeight(select_button_height)
        self.sb_lpno_supersystem_alt.setFixedHeight(select_button_height)

        # File selection buttons for supersystem files, matching the height of the select and header boxes
        sb_lpno_select_supersystem_button = QPushButton("Select Supersystem File")
        sb_lpno_select_supersystem_button.setFixedWidth(150)
        sb_lpno_select_supersystem_button.setFixedHeight(select_button_height)  # Match height to the header

        sb_lpno_select_supersystem_alt_button = QPushButton("Select Supersystem File")
        sb_lpno_select_supersystem_alt_button.setFixedWidth(150)
        sb_lpno_select_supersystem_alt_button.setFixedHeight(select_button_height)  # Match height to the header

        sb_lpno_select_supersystem_button.clicked.connect(lambda: select_file(self.sb_lpno_supersystem_main, "main", 0, "SB_LPNO", is_supersystem=True))
        sb_lpno_select_supersystem_alt_button.clicked.connect(lambda: select_file(self.sb_lpno_supersystem_alt, "alt", 0, "SB_LPNO", is_supersystem=True))

        # Add widgets to the supersystem layout
        sb_lpno_supersystem_layout.addWidget(sb_lpno_select_supersystem_button)
        sb_lpno_supersystem_layout.addWidget(self.sb_lpno_supersystem_main)
        sb_lpno_supersystem_layout.addWidget(sb_lpno_select_supersystem_alt_button)
        sb_lpno_supersystem_layout.addWidget(self.sb_lpno_supersystem_alt)

        layout.addLayout(sb_lpno_supersystem_layout)

        # Handle manual input for supersystem paths in the SB_LPNO section
        self.sb_lpno_supersystem_main.editingFinished.connect(lambda: self.update_filename(self.main_filenames_SB_LPNO, self.normalize_path_gui(self.sb_lpno_supersystem_main.text().strip()), 0))
        self.sb_lpno_supersystem_alt.editingFinished.connect(lambda: self.update_filename(self.alternative_filenames_SB_LPNO, self.normalize_path_gui(self.sb_lpno_supersystem_alt.text().strip()), 0))

        # Subsystem layout for Looser PNO (below supersystem)
        sb_lpno_subsystem_layout = QVBoxLayout()
        layout.addLayout(sb_lpno_subsystem_layout)

        # Button to add subsystem files dynamically
        add_sb_lpno_button = QPushButton("Add Subsystem Files")
        add_sb_lpno_button.setFixedHeight(select_button_height)
        add_sb_lpno_button.clicked.connect(lambda: add_subsystems(sb_lpno_subsystem_layout, "SB_LPNO"))
        layout.addWidget(add_sb_lpno_button)

        layout.addSpacing(10)

        # Adding the "Smaller Basis Set and Tighter PNO" Section
        sb_tpno_label = QLabel("Smaller Basis Set and Tighter PNO")
        sb_tpno_label.setStyleSheet("background-color: #d3e2f4; font-size: 12px; padding: 5px; text-align: center; font-weight: bold;")
        sb_tpno_label.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        layout.addWidget(sb_tpno_label)

        # Add relabel mapping section right after the heading
        self.add_relabel_mapping_section(layout, "SB_TPNO")

        # Create header labels for main and alternative sections
        main_label = QLabel("Main ORCA Output Files")
        alt_label = QLabel("Alternative ORCA Output Files (Optional)")
        main_label.setStyleSheet("background-color: #ebedef; font-size: 11px; padding: 5px; text-align: center; font-weight: bold;")
        alt_label.setStyleSheet("background-color: #ebedef; font-size: 11px; padding: 5px; text-align: center; font-weight: bold;")
        main_label.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        alt_label.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        main_label.setFixedHeight(select_button_height)
        alt_label.setFixedHeight(select_button_height)

        sb_tpno_header_layout = QHBoxLayout()
        sb_tpno_header_layout.addWidget(main_label)
        sb_tpno_header_layout.addSpacerItem(header_spacer)
        sb_tpno_header_layout.addWidget(alt_label)
        layout.addLayout(sb_tpno_header_layout)

        sb_tpno_supersystem_layout = QHBoxLayout()
        self.sb_tpno_supersystem_main = QLineEdit()
        self.sb_tpno_supersystem_alt = QLineEdit()
        self.sb_tpno_supersystem_main.setPlaceholderText("Supersystem Main File")
        self.sb_tpno_supersystem_alt.setPlaceholderText("Supersystem Alternative File")
        self.sb_tpno_supersystem_main.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.sb_tpno_supersystem_alt.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.sb_tpno_supersystem_main.setFixedHeight(select_button_height)
        self.sb_tpno_supersystem_alt.setFixedHeight(select_button_height)

        sb_tpno_select_supersystem_button = QPushButton("Select Supersystem File")
        sb_tpno_select_supersystem_button.setFixedWidth(150)
        sb_tpno_select_supersystem_button.setFixedHeight(select_button_height)

        sb_tpno_select_supersystem_alt_button = QPushButton("Select Supersystem File")
        sb_tpno_select_supersystem_alt_button.setFixedWidth(150)
        sb_tpno_select_supersystem_alt_button.setFixedHeight(select_button_height)

        sb_tpno_select_supersystem_button.clicked.connect(lambda: select_file(self.sb_tpno_supersystem_main, "main", 0, "SB_TPNO", is_supersystem=True))
        sb_tpno_select_supersystem_alt_button.clicked.connect(lambda: select_file(self.sb_tpno_supersystem_alt, "alt", 0, "SB_TPNO", is_supersystem=True))

        sb_tpno_supersystem_layout.addWidget(sb_tpno_select_supersystem_button)
        sb_tpno_supersystem_layout.addWidget(self.sb_tpno_supersystem_main)
        sb_tpno_supersystem_layout.addWidget(sb_tpno_select_supersystem_alt_button)
        sb_tpno_supersystem_layout.addWidget(self.sb_tpno_supersystem_alt)

        layout.addLayout(sb_tpno_supersystem_layout)

        # Handle manual input for the SB_TPNO section
        self.sb_tpno_supersystem_main.editingFinished.connect(lambda: self.update_filename(self.main_filenames_SB_TPNO, self.normalize_path_gui(self.sb_tpno_supersystem_main.text().strip()), 0))
        self.sb_tpno_supersystem_alt.editingFinished.connect(lambda: self.update_filename(self.alternative_filenames_SB_TPNO, self.normalize_path_gui(self.sb_tpno_supersystem_alt.text().strip()), 0))

        sb_tpno_subsystem_layout = QVBoxLayout()
        layout.addLayout(sb_tpno_subsystem_layout)

        # Button to add subsystem files dynamically
        add_sb_tpno_button = QPushButton("Add Subsystem Files")
        add_sb_tpno_button.setFixedHeight(select_button_height)
        add_sb_tpno_button.clicked.connect(lambda: add_subsystems(sb_tpno_subsystem_layout, "SB_TPNO"))
        layout.addWidget(add_sb_tpno_button)

        layout.addSpacing(10)

        # Adding the "Larger Basis Set and Looser PNO" Section
        sb_lpno_label = QLabel("Larger Basis Set and Looser PNO")
        sb_lpno_label.setStyleSheet("background-color: #d3e2f4; font-size: 12px; padding: 5px; text-align: center; font-weight: bold;")
        sb_lpno_label.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        layout.addWidget(sb_lpno_label)

        # Add relabel mapping section right after the heading
        self.add_relabel_mapping_section(layout, "LB_LPNO")

        # Duplicate the header and layout logic for the new section
        main_label = QLabel("Main ORCA Output Files")
        alt_label = QLabel("Alternative ORCA Output Files (Optional)")
        main_label.setStyleSheet("background-color: #ebedef; font-size: 11px; padding: 5px; text-align: center; font-weight: bold;")
        alt_label.setStyleSheet("background-color: #ebedef; font-size: 11px; padding: 5px; text-align: center; font-weight: bold;")
        main_label.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        alt_label.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        main_label.setFixedHeight(select_button_height)
        alt_label.setFixedHeight(select_button_height)

        lb_lpno_header_layout = QHBoxLayout()
        lb_lpno_header_layout.addWidget(main_label)
        lb_lpno_header_layout.addSpacerItem(header_spacer)
        lb_lpno_header_layout.addWidget(alt_label)
        layout.addLayout(lb_lpno_header_layout)

        lb_lpno_supersystem_layout = QHBoxLayout()
        self.lb_lpno_supersystem_main = QLineEdit()
        self.lb_lpno_supersystem_alt = QLineEdit()
        self.lb_lpno_supersystem_main.setPlaceholderText("Supersystem Main File")
        self.lb_lpno_supersystem_alt.setPlaceholderText("Supersystem Alternative File")
        self.lb_lpno_supersystem_main.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.lb_lpno_supersystem_alt.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.lb_lpno_supersystem_main.setFixedHeight(select_button_height)
        self.lb_lpno_supersystem_alt.setFixedHeight(select_button_height)

        lb_lpno_select_supersystem_button = QPushButton("Select Supersystem File")
        lb_lpno_select_supersystem_button.setFixedWidth(150)
        lb_lpno_select_supersystem_button.setFixedHeight(select_button_height)

        lb_lpno_select_supersystem_alt_button = QPushButton("Select Supersystem File")
        lb_lpno_select_supersystem_alt_button.setFixedWidth(150)
        lb_lpno_select_supersystem_alt_button.setFixedHeight(select_button_height)

        lb_lpno_select_supersystem_button.clicked.connect(lambda: select_file(self.lb_lpno_supersystem_main, "main", 0, "LB_LPNO", is_supersystem=True))
        lb_lpno_select_supersystem_alt_button.clicked.connect(lambda: select_file(self.lb_lpno_supersystem_alt, "alt", 0, "LB_LPNO", is_supersystem=True))

        lb_lpno_supersystem_layout.addWidget(lb_lpno_select_supersystem_button)
        lb_lpno_supersystem_layout.addWidget(self.lb_lpno_supersystem_main)
        lb_lpno_supersystem_layout.addWidget(lb_lpno_select_supersystem_alt_button)
        lb_lpno_supersystem_layout.addWidget(self.lb_lpno_supersystem_alt)

        layout.addLayout(lb_lpno_supersystem_layout)

        # Handle manual input for the LB_LPNO section
        self.lb_lpno_supersystem_main.editingFinished.connect(lambda: self.update_filename(self.main_filenames_LB_LPNO, self.normalize_path_gui(self.lb_lpno_supersystem_main.text().strip()), 0))
        self.lb_lpno_supersystem_alt.editingFinished.connect(lambda: self.update_filename(self.alternative_filenames_LB_LPNO, self.normalize_path_gui(self.lb_lpno_supersystem_alt.text().strip()), 0))

        lb_lpno_subsystem_layout = QVBoxLayout()
        layout.addLayout(lb_lpno_subsystem_layout)

        # Button to add subsystem files dynamically
        add_lb_lpno_button = QPushButton("Add Subsystem Files")
        add_lb_lpno_button.setFixedHeight(select_button_height)
        add_lb_lpno_button.clicked.connect(lambda: add_subsystems(lb_lpno_subsystem_layout, "LB_LPNO"))
        layout.addWidget(add_lb_lpno_button)

        layout.addSpacing(10)

        # Adding the "Larger Basis Set and Tighter PNO" Section
        lb_tpno_label = QLabel("Larger Basis Set and Tighter PNO")
        lb_tpno_label.setStyleSheet("background-color: #d3e2f4; font-size: 12px; padding: 5px; text-align: center; font-weight: bold;")
        lb_tpno_label.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        layout.addWidget(lb_tpno_label)

        # Add relabel mapping section right after the heading
        self.add_relabel_mapping_section(layout, "LB_TPNO")

        # Duplicate the header and layout logic for the new section
        main_label = QLabel("Main ORCA Output Files")
        alt_label = QLabel("Alternative ORCA Output Files (Optional)")
        main_label.setStyleSheet("background-color: #ebedef; font-size: 11px; padding: 5px; text-align: center; font-weight: bold;")
        alt_label.setStyleSheet("background-color: #ebedef; font-size: 11px; padding: 5px; text-align: center; font-weight: bold;")
        main_label.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        alt_label.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        main_label.setFixedHeight(select_button_height)
        alt_label.setFixedHeight(select_button_height)

        lb_tpno_header_layout = QHBoxLayout()
        lb_tpno_header_layout.addWidget(main_label)
        lb_tpno_header_layout.addSpacerItem(header_spacer)
        lb_tpno_header_layout.addWidget(alt_label)
        layout.addLayout(lb_tpno_header_layout)

        lb_tpno_supersystem_layout = QHBoxLayout()
        self.lb_tpno_supersystem_main = QLineEdit()
        self.lb_tpno_supersystem_alt = QLineEdit()
        self.lb_tpno_supersystem_main.setPlaceholderText("Supersystem Main File")
        self.lb_tpno_supersystem_alt.setPlaceholderText("Supersystem Alternative File")
        self.lb_tpno_supersystem_main.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.lb_tpno_supersystem_alt.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.lb_tpno_supersystem_main.setFixedHeight(select_button_height)
        self.lb_tpno_supersystem_alt.setFixedHeight(select_button_height)

        lb_tpno_select_supersystem_button = QPushButton("Select Supersystem File")
        lb_tpno_select_supersystem_button.setFixedWidth(150)
        lb_tpno_select_supersystem_button.setFixedHeight(select_button_height)

        lb_tpno_select_supersystem_alt_button = QPushButton("Select Supersystem File")
        lb_tpno_select_supersystem_alt_button.setFixedWidth(150)
        lb_tpno_select_supersystem_alt_button.setFixedHeight(select_button_height)

        lb_tpno_select_supersystem_button.clicked.connect(lambda: select_file(self.lb_tpno_supersystem_main, "main", 0, "LB_TPNO"))
        lb_tpno_select_supersystem_alt_button.clicked.connect(lambda: select_file(self.lb_tpno_supersystem_alt, "alt", 0, "LB_TPNO"))

        lb_tpno_supersystem_layout.addWidget(lb_tpno_select_supersystem_button)
        lb_tpno_supersystem_layout.addWidget(self.lb_tpno_supersystem_main)
        lb_tpno_supersystem_layout.addWidget(lb_tpno_select_supersystem_alt_button)
        lb_tpno_supersystem_layout.addWidget(self.lb_tpno_supersystem_alt)

        layout.addLayout(lb_tpno_supersystem_layout)

        # Handle manual input for the LB_TPNO section
        self.lb_tpno_supersystem_main.editingFinished.connect(lambda: self.update_filename(self.main_filenames_LB_TPNO, self.normalize_path_gui(self.lb_tpno_supersystem_main.text().strip()), 0))
        self.lb_tpno_supersystem_alt.editingFinished.connect(lambda: self.update_filename(self.alternative_filenames_LB_TPNO, self.normalize_path_gui(self.lb_tpno_supersystem_alt.text().strip()), 0))

        lb_tpno_subsystem_layout = QVBoxLayout()
        layout.addLayout(lb_tpno_subsystem_layout)

        # Button to add subsystem files dynamically
        add_lb_tpno_button = QPushButton("Add Subsystem Files")
        add_lb_tpno_button.setFixedHeight(select_button_height)
        add_lb_tpno_button.clicked.connect(lambda: add_subsystems(lb_tpno_subsystem_layout, "LB_TPNO"))
        layout.addWidget(add_lb_tpno_button)

        layout.addSpacing(10)

        # Final spacing and confirm button as before
        layout.addSpacing(10)
        confirm_button = QPushButton("Confirm Files to Proceed")
        confirm_button.setFixedHeight(40)
        confirm_button.setStyleSheet("QPushButton {background-color: #ebedef;} QPushButton:hover {background-color: paleturquoise;}")
        confirm_button.clicked.connect(self.confirm_cps_cbs_nbody_to_proceed)
        layout.addWidget(confirm_button)

        layout.addStretch()

        # Status label for the job process (error message can be copied)
        status_label = QLabel('')
        status_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        status_label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        layout.addWidget(status_label)

        return self.tabs.indexOf(self.nbody_tab)


    def confirm_cps_cbs_nbody_to_proceed(self):
        """Confirm the selected files for all sections, ensuring consistent subsystem counts."""

        # Force an update of the supersystem filenames to capture any last-minute changes from the user input fields
        self.main_filenames_SB_LPNO[0] = self.normalize_path_gui(self.sb_lpno_supersystem_main.text().strip() or '')
        self.alternative_filenames_SB_LPNO[0] = self.normalize_path_gui(self.sb_lpno_supersystem_alt.text().strip() or '')

        self.main_filenames_SB_TPNO[0] = self.normalize_path_gui(self.sb_tpno_supersystem_main.text().strip() or '')
        self.alternative_filenames_SB_TPNO[0] = self.normalize_path_gui(self.sb_tpno_supersystem_alt.text().strip() or '')

        self.main_filenames_LB_LPNO[0] = self.normalize_path_gui(self.lb_lpno_supersystem_main.text().strip() or '')
        self.alternative_filenames_LB_LPNO[0] = self.normalize_path_gui(self.lb_lpno_supersystem_alt.text().strip() or '')

        self.main_filenames_LB_TPNO[0] = self.normalize_path_gui(self.lb_tpno_supersystem_main.text().strip() or '')
        self.alternative_filenames_LB_TPNO[0] = self.normalize_path_gui(self.lb_tpno_supersystem_alt.text().strip() or '')

        # Filter out empty entries from each section and validate without overwriting the lists if validation fails
        validated_main, validated_alt = self.filter_empty_rows_and_validate_main(self.main_filenames_SB_LPNO, self.alternative_filenames_SB_LPNO, "Smaller Basis Set Looser PNO")
        if validated_main is None:  # Return early if validation failed
            return
        self.main_filenames_SB_LPNO, self.alternative_filenames_SB_LPNO = validated_main, validated_alt

        validated_main, validated_alt = self.filter_empty_rows_and_validate_main(self.main_filenames_SB_TPNO, self.alternative_filenames_SB_TPNO, "Smaller Basis Set Tighter PNO")
        if validated_main is None:  # Return early if validation failed
            return
        self.main_filenames_SB_TPNO, self.alternative_filenames_SB_TPNO = validated_main, validated_alt

        validated_main, validated_alt = self.filter_empty_rows_and_validate_main(self.main_filenames_LB_LPNO, self.alternative_filenames_LB_LPNO, "Larger Basis Set Looser PNO")
        if validated_main is None:  # Return early if validation failed
            return
        self.main_filenames_LB_LPNO, self.alternative_filenames_LB_LPNO = validated_main, validated_alt

        validated_main, validated_alt = self.filter_empty_rows_and_validate_main(self.main_filenames_LB_TPNO, self.alternative_filenames_LB_TPNO, "Larger Basis Set Tighter PNO")
        if validated_main is None:  # Return early if validation failed
            return
        self.main_filenames_LB_TPNO, self.alternative_filenames_LB_TPNO = validated_main, validated_alt

        # Count the number of subsystems in each section
        num_subsystems_SB_LPNO = self.count_valid_subsystems(self.main_filenames_SB_LPNO)
        num_subsystems_SB_TPNO = self.count_valid_subsystems(self.main_filenames_SB_TPNO)
        num_subsystems_LB_LPNO = self.count_valid_subsystems(self.main_filenames_LB_LPNO)
        num_subsystems_LB_TPNO = self.count_valid_subsystems(self.main_filenames_LB_TPNO)

        # Check that the number of subsystems is consistent across all sections
        if not (num_subsystems_SB_LPNO == num_subsystems_SB_TPNO == num_subsystems_LB_LPNO == num_subsystems_LB_TPNO):
            QMessageBox.warning(self, "Subsystem Mismatch",
                                "The number of selected systems must be the same across all sections.")
            return  # Stop confirmation if there's a mismatch

        # Ensure at least two valid subsystems are provided in each section
        if not (num_subsystems_SB_LPNO >= 2 and num_subsystems_SB_TPNO >= 2 and num_subsystems_LB_LPNO >= 2 and num_subsystems_LB_TPNO >= 2):
            QMessageBox.warning(self, "Invalid Input", "Please provide one supersystem and at least two subsystem main files in each section.")
            return

        # Prepare the lists of files for confirmation (alternative files can be empty)
        message = "If you want to continue with the following ORCA output files,<br>" \
                  "please press \"OK\". The N-body tab will then be locked to proceed."

        files_overview = f"""<br>
            <b>Smaller Basis Set and Looser PNO</b><br>
            <b>Main Files:</b> {self.main_filenames_SB_LPNO}<br>
            <b>Alternative Files:</b> {self.alternative_filenames_SB_LPNO}<br>
            <b>Relabel Mapping:</b> {self.relabel_mappings.get("SB_LPNO", "Not provided")}<br><br>

            <b>Smaller Basis Set and Tighter PNO</b><br>
            <b>Main Files:</b> {self.main_filenames_SB_TPNO}<br>
            <b>Alternative Files:</b> {self.alternative_filenames_SB_TPNO}<br>
            <b>Relabel Mapping:</b> {self.relabel_mappings.get("SB_TPNO", "Not provided")}<br><br>

            <b>Larger Basis Set and Looser PNO</b><br>
            <b>Main Files:</b> {self.main_filenames_LB_LPNO}<br>
            <b>Alternative Files:</b> {self.alternative_filenames_LB_LPNO}<br>
            <b>Relabel Mapping:</b> {self.relabel_mappings.get("LB_LPNO", "Not provided")}<br><br>

            <b>Larger Basis Set and Tighter PNO</b><br>
            <b>Main Files:</b> {self.main_filenames_LB_TPNO}<br>
            <b>Alternative Files:</b> {self.alternative_filenames_LB_TPNO}<br>
            <b>Relabel Mapping:</b> {self.relabel_mappings.get("LB_TPNO", "Not provided")}<br>
        """

        # Create a scrollable confirmation dialog with selectable text
        confirmation_dialog = QDialog(self)
        confirmation_dialog.setWindowTitle("Confirm Files")

        # Remove the "?" help button from the dialog window
        confirmation_dialog.setWindowFlags(confirmation_dialog.windowFlags() & ~Qt.WindowType.WindowContextHelpButtonHint)

        layout = QVBoxLayout()

        # Add a QLabel for the unselectable portion (the instructional message)
        label = QLabel(message)
        layout.addWidget(label)

        # Add QTextEdit for selectable text and set it to read-only
        text_area = QTextEdit()
        text_area.setReadOnly(True)
        text_area.setHtml(files_overview)  # Use setHtml to apply rich text (bold)

        # Allow the text area to dynamically resize with the window
        text_area.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        layout.addWidget(text_area)

        # Add OK and Cancel buttons
        button_box = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        button_box.accepted.connect(confirmation_dialog.accept)
        button_box.rejected.connect(confirmation_dialog.reject)
        layout.addWidget(button_box)

        confirmation_dialog.setLayout(layout)

        # Adjust the dialog size to fit the contents better, and let it dynamically resize
        confirmation_dialog.resize(670, 400)

        # Show the dialog and check user response
        if confirmation_dialog.exec() == QDialog.DialogCode.Accepted:
            self.lock_tab(self.nbody_tab)
            self.switch_to_next_dynamic_tab()


    def add_cps_nbody_layout(self, layout):
        """Create a CPS layout for the N-body tab"""

        self.nbody_tab.setMinimumSize(700, 400)

        # Initialize lists to hold file paths
        self.main_filenames_LPNO = ['']
        self.alternative_filenames_LPNO = ['']
        self.main_filenames_TPNO = ['']
        self.alternative_filenames_TPNO = ['']

        # Info button to be placed top-right
        info_button = QPushButton('i')
        info_button.setFixedSize(25, 25)
        info_button.setStyleSheet("border-radius: 12px; background-color: lightgray; color: black; font-weight: bold; font-size: 13px;")


        def show_info():
            QMessageBox.information(self, "Alternative ORCA Output File Info",
                                    "If part of the LED data is missing in an ORCA output file but present in another ORCA output file for the same system, both files can be specified as complementary to each other.\n\n"
                                    "For example, if DLPNO-CCSD(T) is chosen, the DLPNO-CCSD output can be used as the main file, while an alternative file that contains the triples contribution but lacks the HF/LED decomposition is specified. Similarly, an HFLD/LED output can be provided as the alternative file for the HF/LED portion.\n\n"
                                    "Note 1: If data for an LED component is present in both the main and alternative files, the main file will be prioritized.\n\n"
                                    "Note 2: As an example, if you perform the calculations with TCutPNO= 1e-6 and TCutPNO=1e-7, then:\n\n"
                                    "Looser PNO: TCutPNO = 1e-6 ; Tighter PNO: TCutPNO = 1e-7")

        info_button.clicked.connect(show_info)

        # Add the info button layout above the header layout
        info_layout = QHBoxLayout()
        info_layout.addStretch()  # Push the info button to the right
        info_layout.addWidget(info_button)
        layout.addLayout(info_layout)  # Info button on its own row

        # Spacer between Main and Alternative headers with a fixed size of 3px
        header_spacer = QSpacerItem(3, 0, QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Minimum)


        def select_file(button, file_type, row, section, is_supersystem=False):
            file_path, _ = QFileDialog.getOpenFileName(self, "Select ORCA Output File", self.last_selected_dir)

            # Check if the user pressed cancel (file_path will be empty)
            if not file_path:
                return  # Do nothing if no file was selected (i.e., cancel was pressed)

            # Normalize the selected file path
            file_path = self.normalize_path_gui(file_path)

            # Update the last selected directory
            self.last_selected_dir = file_path.rsplit("/", 1)[0]

            # Use the ensure_list_size method from the class
            if section == "LPNO":
                self.ensure_list_size(self.main_filenames_LPNO, row)
                self.ensure_list_size(self.alternative_filenames_LPNO, row)
                if file_type == "main":
                    self.main_filenames_LPNO[row] = file_path
                else:
                    self.alternative_filenames_LPNO[row] = file_path or ''

            elif section == "TPNO":
                self.ensure_list_size(self.main_filenames_TPNO, row)
                self.ensure_list_size(self.alternative_filenames_TPNO, row)
                if file_type == "main":
                    self.main_filenames_TPNO[row] = file_path
                else:
                    self.alternative_filenames_TPNO[row] = file_path or ''

            # Update the button text to show the selected file path
            button.setText(file_path)


        def add_subsystems(target_layout, section):
            """Add and process subsystem file paths"""


            def select_multiple_files(section, current_index):
                """Select multiple files and create rows dynamically."""
                file_paths, _ = QFileDialog.getOpenFileNames(self, "Select ORCA Output Files", self.last_selected_dir)

                # Check if the user pressed cancel (file_paths will be empty)
                if not file_paths:
                    return  # Do nothing if no file was selected (i.e., cancel was pressed)

                # Normalize all selected file paths
                file_paths = [self.normalize_path_gui(path) for path in file_paths]

                # Update the last selected directory
                self.last_selected_dir = file_paths[0].rsplit("/", 1)[0] if file_paths else self.last_selected_dir

                # Process each selected file as if the "Add Subsystem Files" button was pressed for each
                for file_path in file_paths:
                    add_single_subsystem_row(file_path, section, current_index)
                    current_index += 1


            def add_single_subsystem_row(file_path, section, row):
                """Create a single row for a subsystem with main and alternative file options."""
                subsystem_row = QHBoxLayout()

                subsystem_main = QLineEdit()
                subsystem_alt = QLineEdit()

                subsystem_main.setPlaceholderText("Subsystem Main File")
                subsystem_alt.setPlaceholderText("Subsystem Alternative File")

                subsystem_main.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
                subsystem_alt.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
                subsystem_main.setFixedHeight(30)
                subsystem_alt.setFixedHeight(30)

                subsystem_main_button = QPushButton("Select Subsystem File")
                subsystem_main_button.setFixedWidth(150)
                subsystem_main_button.setFixedHeight(30)

                subsystem_alt_button = QPushButton("Select Subsystem File")
                subsystem_alt_button.setFixedWidth(150)
                subsystem_alt_button.setFixedHeight(30)

                # Automatically populate the main file with the selected file path
                subsystem_main.setText(file_path)

                if section == "LPNO":
                    self.ensure_list_size(self.main_filenames_LPNO, row)
                    self.ensure_list_size(self.alternative_filenames_LPNO, row)
                    self.main_filenames_LPNO[row] = file_path

                    subsystem_main_button.clicked.connect(lambda: select_file(subsystem_main, "main", row, section))
                    subsystem_alt_button.clicked.connect(lambda: select_file(subsystem_alt, "alt", row, section))

                    subsystem_main.editingFinished.connect(lambda: self.update_filename(self.main_filenames_LPNO, subsystem_main.text(), row))
                    subsystem_alt.editingFinished.connect(lambda: self.update_filename(self.alternative_filenames_LPNO, subsystem_alt.text(), row))

                elif section == "TPNO":
                    self.ensure_list_size(self.main_filenames_TPNO, row)
                    self.ensure_list_size(self.alternative_filenames_TPNO, row)
                    self.main_filenames_TPNO[row] = file_path

                    subsystem_main_button.clicked.connect(lambda: select_file(subsystem_main, "main", row, section))
                    subsystem_alt_button.clicked.connect(lambda: select_file(subsystem_alt, "alt", row, section))

                    subsystem_main.editingFinished.connect(lambda: self.update_filename(self.main_filenames_TPNO, subsystem_main.text(), row))
                    subsystem_alt.editingFinished.connect(lambda: self.update_filename(self.alternative_filenames_TPNO, subsystem_alt.text(), row))

                # Add the buttons and text fields to the layout
                subsystem_row.addWidget(subsystem_main_button)
                subsystem_row.addWidget(subsystem_main)
                subsystem_row.addWidget(subsystem_alt_button)
                subsystem_row.addWidget(subsystem_alt)

                # Add the row to the target layout
                target_layout.addLayout(subsystem_row)

            # Get the current number of subsystems to maintain indexing consistency
            current_index = len(self.main_filenames_LPNO) if section == "LPNO" else len(self.main_filenames_TPNO)

            # Trigger the file selection and row creation process
            select_multiple_files(section, current_index)

        # Adding the "Looser PNO" Section
        lpno_label = QLabel("Looser PNO")
        lpno_label.setStyleSheet("background-color: #d3e2f4; font-size: 12px; padding: 5px; text-align: center; font-weight: bold;")
        lpno_label.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        layout.addWidget(lpno_label)

        # Add relabel mapping section right after the heading
        self.add_relabel_mapping_section(layout, "LPNO")

        # Create header labels for main and alternative sections
        main_label = QLabel("Main ORCA Output Files")
        alt_label = QLabel("Alternative ORCA Output Files (Optional)")

        # Styling the labels: equal size, bold text (font 2px smaller than select button text), and height matching select boxes
        main_label.setStyleSheet("background-color: #ebedef; font-size: 11px; padding: 5px; text-align: center; font-weight: bold;")
        alt_label.setStyleSheet("background-color: #ebedef; font-size: 11px; padding: 5px; text-align: center; font-weight: bold;")

        # Set equal size policy for dynamic resizing of both headers
        main_label.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        alt_label.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)

        # Match the height of the select boxes and headings (assuming 30px for the select buttons)
        select_button_height = 30
        main_label.setFixedHeight(select_button_height)
        alt_label.setFixedHeight(select_button_height)

        # Layout for the headers
        lpno_header_layout = QHBoxLayout()
        lpno_header_layout.addWidget(main_label)  # Main label dynamically resizes
        lpno_header_layout.addSpacerItem(header_spacer)  # Fixed 3px spacer between labels
        lpno_header_layout.addWidget(alt_label)  # Alt label dynamically resizes

        # Add the header layout to the main layout (below the info button)
        layout.addLayout(lpno_header_layout)

        # Supersystem Section (Fixed-size select buttons, expandable populate fields)
        lpno_supersystem_layout = QHBoxLayout()

        self.lpno_supersystem_main = QLineEdit()
        self.lpno_supersystem_alt = QLineEdit()

        # Set placeholders
        self.lpno_supersystem_main.setPlaceholderText("Supersystem Main File")
        self.lpno_supersystem_alt.setPlaceholderText("Supersystem Alternative File")

        # Set size policy and height to match the heading size
        self.lpno_supersystem_main.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.lpno_supersystem_alt.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.lpno_supersystem_main.setFixedHeight(select_button_height)
        self.lpno_supersystem_alt.setFixedHeight(select_button_height)

        # File selection buttons for supersystem files, matching the height of the select and header boxes
        lpno_select_supersystem_button = QPushButton("Select Supersystem File")
        lpno_select_supersystem_button.setFixedWidth(150)
        lpno_select_supersystem_button.setFixedHeight(select_button_height)  # Match height to the header

        lpno_select_supersystem_alt_button = QPushButton("Select Supersystem File")
        lpno_select_supersystem_alt_button.setFixedWidth(150)
        lpno_select_supersystem_alt_button.setFixedHeight(select_button_height)  # Match height to the header

        lpno_select_supersystem_button.clicked.connect(lambda: select_file(self.lpno_supersystem_main, "main", 0, "LPNO", is_supersystem=True))
        lpno_select_supersystem_alt_button.clicked.connect(lambda: select_file(self.lpno_supersystem_alt, "alt", 0, "LPNO", is_supersystem=True))

        # Add widgets to the supersystem layout
        lpno_supersystem_layout.addWidget(lpno_select_supersystem_button)
        lpno_supersystem_layout.addWidget(self.lpno_supersystem_main)
        lpno_supersystem_layout.addWidget(lpno_select_supersystem_alt_button)
        lpno_supersystem_layout.addWidget(self.lpno_supersystem_alt)

        layout.addLayout(lpno_supersystem_layout)

        # Handle manual input for supersystem paths in the LPNO section
        self.lpno_supersystem_main.editingFinished.connect(lambda: self.update_filename(self.main_filenames_LPNO, self.normalize_path_gui(self.lpno_supersystem_main.text().strip()), 0))
        self.lpno_supersystem_alt.editingFinished.connect(lambda: self.update_filename(self.alternative_filenames_LPNO, self.normalize_path_gui(self.lpno_supersystem_alt.text().strip()), 0))

        # Subsystem layout for Looser PNO (below supersystem)
        lpno_subsystem_layout = QVBoxLayout()
        layout.addLayout(lpno_subsystem_layout)

        # Button to add subsystem files dynamically
        add_lpno_button = QPushButton("Add Subsystem Files")
        add_lpno_button.setFixedHeight(select_button_height)
        add_lpno_button.clicked.connect(lambda: add_subsystems(lpno_subsystem_layout, "LPNO"))
        layout.addWidget(add_lpno_button)

        layout.addSpacing(10)

        # Adding the "Tighter PNO" Section
        tpno_label = QLabel("Tighter PNO")
        tpno_label.setStyleSheet("background-color: #d3e2f4; font-size: 12px; padding: 5px; text-align: center; font-weight: bold;")
        tpno_label.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        layout.addWidget(tpno_label)

        # Add relabel mapping section right after the heading
        self.add_relabel_mapping_section(layout, "TPNO")

        # Create header labels for main and alternative sections
        main_label = QLabel("Main ORCA Output Files")
        alt_label = QLabel("Alternative ORCA Output Files (Optional)")
        main_label.setStyleSheet("background-color: #ebedef; font-size: 11px; padding: 5px; text-align: center; font-weight: bold;")
        alt_label.setStyleSheet("background-color: #ebedef; font-size: 11px; padding: 5px; text-align: center; font-weight: bold;")
        main_label.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        alt_label.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        main_label.setFixedHeight(select_button_height)
        alt_label.setFixedHeight(select_button_height)

        tpno_header_layout = QHBoxLayout()
        tpno_header_layout.addWidget(main_label)
        tpno_header_layout.addSpacerItem(header_spacer)
        tpno_header_layout.addWidget(alt_label)
        layout.addLayout(tpno_header_layout)

        tpno_supersystem_layout = QHBoxLayout()
        self.tpno_supersystem_main = QLineEdit()
        self.tpno_supersystem_alt = QLineEdit()
        self.tpno_supersystem_main.setPlaceholderText("Supersystem Main File")
        self.tpno_supersystem_alt.setPlaceholderText("Supersystem Alternative File")
        self.tpno_supersystem_main.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.tpno_supersystem_alt.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.tpno_supersystem_main.setFixedHeight(select_button_height)
        self.tpno_supersystem_alt.setFixedHeight(select_button_height)

        tpno_select_supersystem_button = QPushButton("Select Supersystem File")
        tpno_select_supersystem_button.setFixedWidth(150)
        tpno_select_supersystem_button.setFixedHeight(select_button_height)

        tpno_select_supersystem_alt_button = QPushButton("Select Supersystem File")
        tpno_select_supersystem_alt_button.setFixedWidth(150)
        tpno_select_supersystem_alt_button.setFixedHeight(select_button_height)

        tpno_select_supersystem_button.clicked.connect(lambda: select_file(self.tpno_supersystem_main, "main", 0, "TPNO", is_supersystem=True))
        tpno_select_supersystem_alt_button.clicked.connect(lambda: select_file(self.tpno_supersystem_alt, "alt", 0, "TPNO", is_supersystem=True))

        tpno_supersystem_layout.addWidget(tpno_select_supersystem_button)
        tpno_supersystem_layout.addWidget(self.tpno_supersystem_main)
        tpno_supersystem_layout.addWidget(tpno_select_supersystem_alt_button)
        tpno_supersystem_layout.addWidget(self.tpno_supersystem_alt)

        layout.addLayout(tpno_supersystem_layout)

        # Handle manual input for the TPNO section
        self.tpno_supersystem_main.editingFinished.connect(lambda: self.update_filename(self.main_filenames_TPNO, self.normalize_path_gui(self.tpno_supersystem_main.text().strip()), 0))
        self.tpno_supersystem_alt.editingFinished.connect(lambda: self.update_filename(self.alternative_filenames_TPNO, self.normalize_path_gui(self.tpno_supersystem_alt.text().strip()), 0))

        tpno_subsystem_layout = QVBoxLayout()
        layout.addLayout(tpno_subsystem_layout)

        # Button to add subsystem files dynamically
        add_tpno_button = QPushButton("Add Subsystem Files")
        add_tpno_button.setFixedHeight(select_button_height)
        add_tpno_button.clicked.connect(lambda: add_subsystems(tpno_subsystem_layout, "TPNO"))
        layout.addWidget(add_tpno_button)

        layout.addSpacing(10)

        # Final spacing and confirm button as before
        layout.addSpacing(10)
        confirm_button = QPushButton("Confirm Files to Proceed")
        confirm_button.setFixedHeight(40)
        confirm_button.setStyleSheet("QPushButton {background-color: #ebedef;} QPushButton:hover {background-color: paleturquoise;}")
        confirm_button.clicked.connect(self.confirm_cps_nbody_to_proceed)
        layout.addWidget(confirm_button)

        layout.addSpacing(120)

        # Status label for the job process (error message can be copied)
        status_label = QLabel('')
        status_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        status_label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        layout.addWidget(status_label)

        return self.tabs.indexOf(self.nbody_tab)


    def confirm_cps_nbody_to_proceed(self):
        """Confirm the selected files for all sections, ensuring consistent subsystem counts."""

        # Force an update of the supersystem filenames to capture any last-minute changes from the user input fields
        self.main_filenames_LPNO[0] = self.normalize_path_gui(self.lpno_supersystem_main.text().strip() or '')
        self.alternative_filenames_LPNO[0] = self.normalize_path_gui(self.lpno_supersystem_alt.text().strip() or '')

        self.main_filenames_TPNO[0] = self.normalize_path_gui(self.tpno_supersystem_main.text().strip() or '')
        self.alternative_filenames_TPNO[0] = self.normalize_path_gui(self.tpno_supersystem_alt.text().strip() or '')

        # Filter out empty entries from each section and validate without overwriting the lists if validation fails
        validated_main, validated_alt = self.filter_empty_rows_and_validate_main(self.main_filenames_LPNO, self.alternative_filenames_LPNO, "Looser PNO")
        if validated_main is None:  # Return early if validation failed
            return
        self.main_filenames_LPNO, self.alternative_filenames_LPNO = validated_main, validated_alt

        validated_main, validated_alt = self.filter_empty_rows_and_validate_main(self.main_filenames_TPNO, self.alternative_filenames_TPNO, "Tighter PNO")
        if validated_main is None:  # Return early if validation failed
            return
        self.main_filenames_TPNO, self.alternative_filenames_TPNO = validated_main, validated_alt

        # Count the number of subsystems in each section
        num_subsystems_LPNO = self.count_valid_subsystems(self.main_filenames_LPNO)
        num_subsystems_TPNO = self.count_valid_subsystems(self.main_filenames_TPNO)

        # Check that the number of subsystems is consistent across all sections
        if not (num_subsystems_LPNO == num_subsystems_TPNO):
            QMessageBox.warning(self, "Subsystem Mismatch",
                                "The number of selected systems must be the same across all sections.")
            return  # Stop confirmation if there's a mismatch

        # Ensure at least two valid subsystems are provided in each section
        if not (num_subsystems_LPNO >= 2 and num_subsystems_TPNO >= 2):
            QMessageBox.warning(self, "Invalid Input", "Please provide one supersystem and at least two subsystem main files in each section.")
            return

        # Prepare the lists of files for confirmation (alternative files can be empty)
        message = "If you want to continue with the following ORCA output files,<br>" \
                  "please press \"OK\". The N-body tab will then be locked to proceed."

        files_overview = f"""<br>
            <b>Looser PNO</b><br>
            <b>Main Files:</b> {self.main_filenames_LPNO}<br>
            <b>Alternative Files:</b> {self.alternative_filenames_LPNO}<br>
            <b>Relabel Mapping:</b> {self.relabel_mappings.get("LPNO", "Not provided")}<br><br>

            <b>Tighter PNO</b><br>
            <b>Main Files:</b> {self.main_filenames_TPNO}<br>
            <b>Alternative Files:</b> {self.alternative_filenames_TPNO}<br>
            <b>Relabel Mapping:</b> {self.relabel_mappings.get("TPNO", "Not provided")}<br>"""

        # Create a scrollable confirmation dialog with selectable text
        confirmation_dialog = QDialog(self)
        confirmation_dialog.setWindowTitle("Confirm Files")

        # Remove the "?" help button from the dialog window
        confirmation_dialog.setWindowFlags(confirmation_dialog.windowFlags() & ~Qt.WindowType.WindowContextHelpButtonHint)

        layout = QVBoxLayout()

        # Add a QLabel for the unselectable portion (the instructional message)
        label = QLabel(message)
        layout.addWidget(label)

        # Add QTextEdit for selectable text and set it to read-only
        text_area = QTextEdit()
        text_area.setReadOnly(True)
        text_area.setHtml(files_overview)  # Use setHtml to apply rich text (bold)

        # Allow the text area to dynamically resize with the window
        text_area.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        layout.addWidget(text_area)

        # Add OK and Cancel buttons
        button_box = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        button_box.accepted.connect(confirmation_dialog.accept)
        button_box.rejected.connect(confirmation_dialog.reject)
        layout.addWidget(button_box)

        confirmation_dialog.setLayout(layout)

        # Adjust the dialog size to fit the contents better, and let it dynamically resize
        confirmation_dialog.resize(670, 400)

        # Show the dialog and check user response
        if confirmation_dialog.exec() == QDialog.DialogCode.Accepted:
            self.lock_tab(self.nbody_tab)
            self.switch_to_next_dynamic_tab()


    def add_cbs_nbody_layout(self, layout):
        """Create a CBS layout for the N-body tab"""

        self.nbody_tab.setMinimumSize(700, 400)

        # Initialize lists to hold file paths
        self.main_filenames_SB = ['']
        self.alternative_filenames_SB = ['']
        self.main_filenames_LB = ['']
        self.alternative_filenames_LB = ['']

        # Info button to be placed top-right
        info_button = QPushButton('i')
        info_button.setFixedSize(25, 25)
        info_button.setStyleSheet("border-radius: 12px; background-color: lightgray; color: black; font-weight: bold; font-size: 13px;")


        def show_info():
            QMessageBox.information(self, "Alternative ORCA Output File Info",
                                    "If part of the LED data is missing in an ORCA output file but present in another ORCA output file for the same system, both files can be specified as complementary to each other.\n\n"
                                    "For example, if DLPNO-CCSD(T) is chosen, the DLPNO-CCSD output can be used as the main file, while an alternative file that contains the triples contribution but lacks the HF/LED decomposition is specified. Similarly, an HFLD/LED output can be provided as the alternative file for the HF/LED portion.\n\n"
                                    "Note 1: If data for an LED component is present in both the main and alternative files, the main file will be prioritized.\n\n"
                                    "Note 2: As an example, if you perform the calculations with cc-pVTZ and cc-pVQZ basis sets, then:\n\n"
                                    "Smaller Basis Set: cc-pVTZ; Larger Basis Set: cc-pVQZ")

        info_button.clicked.connect(show_info)

        # Add the info button layout above the header layout
        info_layout = QHBoxLayout()
        info_layout.addStretch()  # Push the info button to the right
        info_layout.addWidget(info_button)
        layout.addLayout(info_layout)  # Info button on its own row

        # Spacer between Main and Alternative headers with a fixed size of 3px
        header_spacer = QSpacerItem(3, 0, QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Minimum)


        def select_file(button, file_type, row, section, is_supersystem=False):
            file_path, _ = QFileDialog.getOpenFileName(self, "Select ORCA Output File", self.last_selected_dir)

            # Check if the user pressed cancel (file_path will be empty)
            if not file_path:
                return  # Do nothing if no file was selected (i.e., cancel was pressed)

            # Normalize the selected file path
            file_path = self.normalize_path_gui(file_path)

            # Update the last selected directory
            self.last_selected_dir = file_path.rsplit("/", 1)[0]

            # Use the ensure_list_size method from the class
            if section == "SB":
                self.ensure_list_size(self.main_filenames_SB, row)
                self.ensure_list_size(self.alternative_filenames_SB, row)
                if file_type == "main":
                    self.main_filenames_SB[row] = file_path
                else:
                    self.alternative_filenames_SB[row] = file_path or ''

            elif section == "LB":
                self.ensure_list_size(self.main_filenames_LB, row)
                self.ensure_list_size(self.alternative_filenames_LB, row)
                if file_type == "main":
                    self.main_filenames_LB[row] = file_path
                else:
                    self.alternative_filenames_LB[row] = file_path or ''

            # Update the button text to show the selected file path
            button.setText(file_path)


        def add_subsystems(target_layout, section):
            """Add and process subsystem file paths"""


            def select_multiple_files(section, current_index):
                """Select multiple files and create rows dynamically."""
                file_paths, _ = QFileDialog.getOpenFileNames(self, "Select ORCA Output Files", self.last_selected_dir)

                # Check if the user pressed cancel (file_paths will be empty)
                if not file_paths:
                    return  # Do nothing if no file was selected (i.e., cancel was pressed)

                # Normalize all selected file paths
                file_paths = [self.normalize_path_gui(path) for path in file_paths]

                # Update the last selected directory
                self.last_selected_dir = file_paths[0].rsplit("/", 1)[0] if file_paths else self.last_selected_dir

                # Process each selected file as if the "Add Subsystem Files" button was pressed for each
                for file_path in file_paths:
                    add_single_subsystem_row(file_path, section, current_index)
                    current_index += 1


            def add_single_subsystem_row(file_path, section, row):
                """Create a single row for a subsystem with main and alternative file options."""
                subsystem_row = QHBoxLayout()

                subsystem_main = QLineEdit()
                subsystem_alt = QLineEdit()

                subsystem_main.setPlaceholderText("Subsystem Main File")
                subsystem_alt.setPlaceholderText("Subsystem Alternative File")

                subsystem_main.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
                subsystem_alt.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
                subsystem_main.setFixedHeight(30)
                subsystem_alt.setFixedHeight(30)

                subsystem_main_button = QPushButton("Select Subsystem File")
                subsystem_main_button.setFixedWidth(150)
                subsystem_main_button.setFixedHeight(30)

                subsystem_alt_button = QPushButton("Select Subsystem File")
                subsystem_alt_button.setFixedWidth(150)
                subsystem_alt_button.setFixedHeight(30)

                # Automatically populate the main file with the selected file path
                subsystem_main.setText(file_path)

                if section == "SB":
                    self.ensure_list_size(self.main_filenames_SB, row)
                    self.ensure_list_size(self.alternative_filenames_SB, row)
                    self.main_filenames_SB[row] = file_path

                    subsystem_main_button.clicked.connect(lambda: select_file(subsystem_main, "main", row, section))
                    subsystem_alt_button.clicked.connect(lambda: select_file(subsystem_alt, "alt", row, section))

                    subsystem_main.editingFinished.connect(lambda: self.update_filename(self.main_filenames_SB, subsystem_main.text(), row))
                    subsystem_alt.editingFinished.connect(lambda: self.update_filename(self.alternative_filenames_SB, subsystem_alt.text(), row))

                elif section == "LB":
                    self.ensure_list_size(self.main_filenames_LB, row)
                    self.ensure_list_size(self.alternative_filenames_LB, row)
                    self.main_filenames_LB[row] = file_path

                    subsystem_main_button.clicked.connect(lambda: select_file(subsystem_main, "main", row, section))
                    subsystem_alt_button.clicked.connect(lambda: select_file(subsystem_alt, "alt", row, section))

                    subsystem_main.editingFinished.connect(lambda: self.update_filename(self.main_filenames_LB, subsystem_main.text(), row))
                    subsystem_alt.editingFinished.connect(lambda: self.update_filename(self.alternative_filenames_LB, subsystem_alt.text(), row))

                # Add the buttons and text fields to the layout
                subsystem_row.addWidget(subsystem_main_button)
                subsystem_row.addWidget(subsystem_main)
                subsystem_row.addWidget(subsystem_alt_button)
                subsystem_row.addWidget(subsystem_alt)

                # Add the row to the target layout
                target_layout.addLayout(subsystem_row)

            # Get the current number of subsystems to maintain indexing consistency
            current_index = len(self.main_filenames_SB) if section == "SB" else len(self.main_filenames_LB)

            # Trigger the file selection and row creation process
            select_multiple_files(section, current_index)

        # Adding the "Smaller Basis Set" Section
        sb_label = QLabel("Smaller Basis Set")
        sb_label.setStyleSheet("background-color: #d3e2f4; font-size: 12px; padding: 5px; text-align: center; font-weight: bold;")
        sb_label.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        layout.addWidget(sb_label)

        # Add relabel mapping section right after the heading
        self.add_relabel_mapping_section(layout, "SB")

        # Create header labels for main and alternative sections
        main_label = QLabel("Main ORCA Output Files")
        alt_label = QLabel("Alternative ORCA Output Files (Optional)")

        # Styling the labels: equal size, bold text (font 2px smaller than select button text), and height matching select boxes
        main_label.setStyleSheet("background-color: #ebedef; font-size: 11px; padding: 5px; text-align: center; font-weight: bold;")
        alt_label.setStyleSheet("background-color: #ebedef; font-size: 11px; padding: 5px; text-align: center; font-weight: bold;")

        # Set equal size policy for dynamic resizing of both headers
        main_label.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        alt_label.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)

        # Match the height of the select boxes and headings (assuming 30px for the select buttons)
        select_button_height = 30
        main_label.setFixedHeight(select_button_height)
        alt_label.setFixedHeight(select_button_height)

        # Layout for the headers
        sb_header_layout = QHBoxLayout()
        sb_header_layout.addWidget(main_label)  # Main label dynamically resizes
        sb_header_layout.addSpacerItem(header_spacer)  # Fixed 3px spacer between labels
        sb_header_layout.addWidget(alt_label)  # Alt label dynamically resizes

        # Add the header layout to the main layout (below the info button)
        layout.addLayout(sb_header_layout)

        # Supersystem Section (Fixed-size select buttons, expandable populate fields)
        sb_supersystem_layout = QHBoxLayout()

        self.sb_supersystem_main = QLineEdit()
        self.sb_supersystem_alt = QLineEdit()

        # Set placeholders
        self.sb_supersystem_main.setPlaceholderText("Supersystem Main File")
        self.sb_supersystem_alt.setPlaceholderText("Supersystem Alternative File")

        # Set size policy and height to match the heading size
        self.sb_supersystem_main.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.sb_supersystem_alt.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.sb_supersystem_main.setFixedHeight(select_button_height)
        self.sb_supersystem_alt.setFixedHeight(select_button_height)

        # File selection buttons for supersystem files, matching the height of the select and header boxes
        sb_select_supersystem_button = QPushButton("Select Supersystem File")
        sb_select_supersystem_button.setFixedWidth(150)
        sb_select_supersystem_button.setFixedHeight(select_button_height)  # Match height to the header

        sb_select_supersystem_alt_button = QPushButton("Select Supersystem File")
        sb_select_supersystem_alt_button.setFixedWidth(150)
        sb_select_supersystem_alt_button.setFixedHeight(select_button_height)  # Match height to the header

        sb_select_supersystem_button.clicked.connect(lambda: select_file(self.sb_supersystem_main, "main", 0, "SB", is_supersystem=True))
        sb_select_supersystem_alt_button.clicked.connect(lambda: select_file(self.sb_supersystem_alt, "alt", 0, "SB", is_supersystem=True))

        # Add widgets to the supersystem layout
        sb_supersystem_layout.addWidget(sb_select_supersystem_button)
        sb_supersystem_layout.addWidget(self.sb_supersystem_main)
        sb_supersystem_layout.addWidget(sb_select_supersystem_alt_button)
        sb_supersystem_layout.addWidget(self.sb_supersystem_alt)

        layout.addLayout(sb_supersystem_layout)

        # Handle manual input for supersystem paths in the SB section
        self.sb_supersystem_main.editingFinished.connect(lambda: self.update_filename(self.main_filenames_SB, self.normalize_path_gui(self.sb_supersystem_main.text().strip()), 0))
        self.sb_supersystem_alt.editingFinished.connect(lambda: self.update_filename(self.alternative_filenames_SB, self.normalize_path_gui(self.sb_supersystem_alt.text().strip()), 0))

        # Subsystem layout for Smaller Basis Set (below supersystem)
        sb_subsystem_layout = QVBoxLayout()
        layout.addLayout(sb_subsystem_layout)

        # Button to add subsystem files dynamically
        add_sb_button = QPushButton("Add Subsystem Files")
        add_sb_button.setFixedHeight(select_button_height)
        add_sb_button.clicked.connect(lambda: add_subsystems(sb_subsystem_layout, "SB"))
        layout.addWidget(add_sb_button)

        layout.addSpacing(10)

        # Adding the "Larger Basis Set" Section
        lb_label = QLabel("Larger Basis Set")
        lb_label.setStyleSheet("background-color: #d3e2f4; font-size: 12px; padding: 5px; text-align: center; font-weight: bold;")
        lb_label.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        layout.addWidget(lb_label)

        # Add relabel mapping section right after the heading
        self.add_relabel_mapping_section(layout, "LB")

        # Create header labels for main and alternative sections
        main_label = QLabel("Main ORCA Output Files")
        alt_label = QLabel("Alternative ORCA Output Files (Optional)")
        main_label.setStyleSheet("background-color: #ebedef; font-size: 11px; padding: 5px; text-align: center; font-weight: bold;")
        alt_label.setStyleSheet("background-color: #ebedef; font-size: 11px; padding: 5px; text-align: center; font-weight: bold;")
        main_label.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        alt_label.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        main_label.setFixedHeight(select_button_height)
        alt_label.setFixedHeight(select_button_height)

        lb_header_layout = QHBoxLayout()
        lb_header_layout.addWidget(main_label)
        lb_header_layout.addSpacerItem(header_spacer)
        lb_header_layout.addWidget(alt_label)
        layout.addLayout(lb_header_layout)

        lb_supersystem_layout = QHBoxLayout()
        self.lb_supersystem_main = QLineEdit()
        self.lb_supersystem_alt = QLineEdit()
        self.lb_supersystem_main.setPlaceholderText("Supersystem Main File")
        self.lb_supersystem_alt.setPlaceholderText("Supersystem Alternative File")
        self.lb_supersystem_main.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.lb_supersystem_alt.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.lb_supersystem_main.setFixedHeight(select_button_height)
        self.lb_supersystem_alt.setFixedHeight(select_button_height)

        lb_select_supersystem_button = QPushButton("Select Supersystem File")
        lb_select_supersystem_button.setFixedWidth(150)
        lb_select_supersystem_button.setFixedHeight(select_button_height)

        lb_select_supersystem_alt_button = QPushButton("Select Supersystem File")
        lb_select_supersystem_alt_button.setFixedWidth(150)
        lb_select_supersystem_alt_button.setFixedHeight(select_button_height)

        lb_select_supersystem_button.clicked.connect(lambda: select_file(self.lb_supersystem_main, "main", 0, "LB", is_supersystem=True))
        lb_select_supersystem_alt_button.clicked.connect(lambda: select_file(self.lb_supersystem_alt, "alt", 0, "LB", is_supersystem=True))

        lb_supersystem_layout.addWidget(lb_select_supersystem_button)
        lb_supersystem_layout.addWidget(self.lb_supersystem_main)
        lb_supersystem_layout.addWidget(lb_select_supersystem_alt_button)
        lb_supersystem_layout.addWidget(self.lb_supersystem_alt)

        layout.addLayout(lb_supersystem_layout)

        # Handle manual input for the LB section
        self.lb_supersystem_main.editingFinished.connect(lambda: self.update_filename(self.main_filenames_LB, self.normalize_path_gui(self.lb_supersystem_main.text().strip()), 0))
        self.lb_supersystem_alt.editingFinished.connect(lambda: self.update_filename(self.alternative_filenames_LB, self.normalize_path_gui(self.lb_supersystem_alt.text().strip()), 0))

        lb_subsystem_layout = QVBoxLayout()
        layout.addLayout(lb_subsystem_layout)

        # Button to add subsystem files dynamically
        add_lb_button = QPushButton("Add Subsystem Files")
        add_lb_button.setFixedHeight(select_button_height)
        add_lb_button.clicked.connect(lambda: add_subsystems(lb_subsystem_layout, "LB"))
        layout.addWidget(add_lb_button)

        layout.addSpacing(10)

        # Final spacing and confirm button as before
        layout.addSpacing(10)
        confirm_button = QPushButton("Confirm Files to Proceed")
        confirm_button.setFixedHeight(40)
        confirm_button.setStyleSheet("QPushButton {background-color: #ebedef;} QPushButton:hover {background-color: paleturquoise;}")
        confirm_button.clicked.connect(self.confirm_cbs_nbody_to_proceed)
        layout.addWidget(confirm_button)

        layout.addSpacing(120)

        # Status label for the job process (error message can be copied)
        status_label = QLabel('')
        status_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        status_label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        layout.addWidget(status_label)

        return self.tabs.indexOf(self.nbody_tab)


    def confirm_cbs_nbody_to_proceed(self):
        """Confirm the selected files for all sections, ensuring consistent subsystem counts."""

        # Force an update of the supersystem filenames to capture any last-minute changes from the user input fields
        self.main_filenames_SB[0] = self.normalize_path_gui(self.sb_supersystem_main.text().strip() or '')
        self.alternative_filenames_SB[0] = self.normalize_path_gui(self.sb_supersystem_alt.text().strip() or '')

        self.main_filenames_LB[0] = self.normalize_path_gui(self.lb_supersystem_main.text().strip() or '')
        self.alternative_filenames_LB[0] = self.normalize_path_gui(self.lb_supersystem_alt.text().strip() or '')

        # Filter out empty entries from each section and validate without overwriting the lists if validation fails
        validated_main, validated_alt = self.filter_empty_rows_and_validate_main(self.main_filenames_SB, self.alternative_filenames_SB, "Smaller Basis Set")
        if validated_main is None:  # Return early if validation failed
            return
        self.main_filenames_SB, self.alternative_filenames_SB = validated_main, validated_alt

        validated_main, validated_alt = self.filter_empty_rows_and_validate_main(self.main_filenames_LB, self.alternative_filenames_LB, "Larger Basis Set")
        if validated_main is None:  # Return early if validation failed
            return
        self.main_filenames_LB, self.alternative_filenames_LB = validated_main, validated_alt

        # Count the number of subsystems in each section
        num_subsystems_SB = self.count_valid_subsystems(self.main_filenames_SB)
        num_subsystems_LB = self.count_valid_subsystems(self.main_filenames_LB)

        # Check that the number of subsystems is consistent across all sections
        if not (num_subsystems_SB == num_subsystems_LB):
            QMessageBox.warning(self, "Subsystem Mismatch",
                                "The number of selected systems must be the same across all sections.")
            return  # Stop confirmation if there's a mismatch

        # Ensure at least two valid subsystems are provided in each section
        if not (num_subsystems_SB >= 2 and num_subsystems_LB >= 2):
            QMessageBox.warning(self, "Invalid Input", "Please provide one supersystem and at least two subsystem main files in each section.")
            return

        # Prepare the lists of files for confirmation (alternative files can be empty)
        message = "If you want to continue with the following ORCA output files,<br>" \
                  "please press \"OK\". The N-body tab will then be locked to proceed."

        files_overview = f"""<br>
            <b>Smaller Basis Set</b><br>
            <b>Main Files:</b> {self.main_filenames_SB}<br>
            <b>Alternative Files:</b> {self.alternative_filenames_SB}<br>
            <b>Relabel Mapping:</b> {self.relabel_mappings.get("SB", "Not provided")}<br><br>

            <b>Larger Basis Set</b><br>
            <b>Main Files:</b> {self.main_filenames_LB}<br>
            <b>Alternative Files:</b> {self.alternative_filenames_LB}<br>
            <b>Relabel Mapping:</b> {self.relabel_mappings.get("LB", "Not provided")}<br>"""

        # Create a scrollable confirmation dialog with selectable text
        confirmation_dialog = QDialog(self)
        confirmation_dialog.setWindowTitle("Confirm Files")

        # Remove the "?" help button from the dialog window
        confirmation_dialog.setWindowFlags(confirmation_dialog.windowFlags() & ~Qt.WindowType.WindowContextHelpButtonHint)

        layout = QVBoxLayout()

        # Add a QLabel for the unselectable portion (the instructional message)
        label = QLabel(message)
        layout.addWidget(label)

        # Add QTextEdit for selectable text and set it to read-only
        text_area = QTextEdit()
        text_area.setReadOnly(True)
        text_area.setHtml(files_overview)  # Use setHtml to apply rich text (bold)

        # Allow the text area to dynamically resize with the window
        text_area.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        layout.addWidget(text_area)

        # Add OK and Cancel buttons
        button_box = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        button_box.accepted.connect(confirmation_dialog.accept)
        button_box.rejected.connect(confirmation_dialog.reject)
        layout.addWidget(button_box)

        confirmation_dialog.setLayout(layout)

        # Adjust the dialog size to fit the contents better, and let it dynamically resize
        confirmation_dialog.resize(670, 400)

        # Show the dialog and check user response
        if confirmation_dialog.exec() == QDialog.DialogCode.Accepted:
            self.lock_tab(self.nbody_tab)
            self.switch_to_next_dynamic_tab()


    def add_standard_nbody_layout(self, layout):
        """Create a standard layout for the N-body tab"""

        self.nbody_tab.setMinimumSize(700, 400)

        # Initialize lists to hold file paths
        self.main_filenames = ['']
        self.alternative_filenames = ['']

        # Info button to be placed top-right
        info_button = QPushButton('i')
        info_button.setFixedSize(25, 25)
        info_button.setStyleSheet("border-radius: 12px; background-color: lightgray; color: black; font-weight: bold; font-size: 13px;")


        def show_info():
            QMessageBox.information(self, "Alternative ORCA Output File Info",
                                    "If part of the LED data is missing in an ORCA output file but present in another ORCA output file for the same system, both files can be specified as complementary to each other.\n\n"
                                    "For example, if DLPNO-CCSD(T) is chosen, the DLPNO-CCSD output can be used as the main file, while an alternative file that contains the triples contribution but lacks the HF/LED decomposition is specified. Similarly, an HFLD/LED output can be provided as the alternative file for the HF/LED portion.\n\n"
                                    "Note: If data for an LED component is present in both the main and alternative files, the main file will be prioritized.")

        info_button.clicked.connect(show_info)

        # Add the info button layout above the header layout
        info_layout = QHBoxLayout()
        info_layout.addStretch()  # Push the info button to the right
        info_layout.addWidget(info_button)
        layout.addLayout(info_layout)  # Info button on its own row

        # Spacer between Main and Alternative headers with a fixed size of 3px
        header_spacer = QSpacerItem(3, 0, QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Minimum)


        def select_file(button, file_type, row, section, is_supersystem=False):
            file_path, _ = QFileDialog.getOpenFileName(self, "Select ORCA Output File", self.last_selected_dir)

            # Check if the user pressed cancel (file_path will be empty)
            if not file_path:
                return  # Do nothing if no file was selected (i.e., cancel was pressed)

            # Normalize the selected file path
            file_path = self.normalize_path_gui(file_path)

            # Update the last selected directory
            self.last_selected_dir = file_path.rsplit("/", 1)[0]

            # Use the ensure_list_size method from the class
            self.ensure_list_size(self.main_filenames, row)
            self.ensure_list_size(self.alternative_filenames, row)
            if file_type == "main":
                self.main_filenames[row] = file_path
            else:
                self.alternative_filenames[row] = file_path or ''


            # Update the button text to show the selected file path
            button.setText(file_path)


        def add_subsystems(target_layout, section):
            """Add and process subsystem file paths"""

            def select_multiple_files(section, current_index):
                """Select multiple files and create rows dynamically."""
                file_paths, _ = QFileDialog.getOpenFileNames(self, "Select ORCA Output Files", self.last_selected_dir)

                # Check if the user pressed cancel (file_paths will be empty)
                if not file_paths:
                    return  # Do nothing if no file was selected (i.e., cancel was pressed)

                # Normalize all selected file paths
                file_paths = [self.normalize_path_gui(path) for path in file_paths]

                # Update the last selected directory
                self.last_selected_dir = file_paths[0].rsplit("/", 1)[0] if file_paths else self.last_selected_dir

                # Process each selected file as if the "Add Subsystem Files" button was pressed for each
                for file_path in file_paths:
                    add_single_subsystem_row(file_path, section, current_index)
                    current_index += 1


            def add_single_subsystem_row(file_path, section, row):
                """Create a single row for a subsystem with main and alternative file options."""
                subsystem_row = QHBoxLayout()

                subsystem_main = QLineEdit()
                subsystem_alt = QLineEdit()

                subsystem_main.setPlaceholderText("Subsystem Main File")
                subsystem_alt.setPlaceholderText("Subsystem Alternative File")

                subsystem_main.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
                subsystem_alt.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
                subsystem_main.setFixedHeight(30)
                subsystem_alt.setFixedHeight(30)

                subsystem_main_button = QPushButton("Select Subsystem File")
                subsystem_main_button.setFixedWidth(150)
                subsystem_main_button.setFixedHeight(30)

                subsystem_alt_button = QPushButton("Select Subsystem File")
                subsystem_alt_button.setFixedWidth(150)
                subsystem_alt_button.setFixedHeight(30)

                # Automatically populate the main file with the selected file path
                subsystem_main.setText(file_path)

                self.ensure_list_size(self.main_filenames, row)
                self.ensure_list_size(self.alternative_filenames, row)
                self.main_filenames[row] = file_path

                subsystem_main_button.clicked.connect(lambda: select_file(subsystem_main, "main", row, section))
                subsystem_alt_button.clicked.connect(lambda: select_file(subsystem_alt, "alt", row, section))

                subsystem_main.editingFinished.connect(lambda: self.update_filename(self.main_filenames, subsystem_main.text(), row))
                subsystem_alt.editingFinished.connect(lambda: self.update_filename(self.alternative_filenames, subsystem_alt.text(), row))

                # Add the buttons and text fields to the layout
                subsystem_row.addWidget(subsystem_main_button)
                subsystem_row.addWidget(subsystem_main)
                subsystem_row.addWidget(subsystem_alt_button)
                subsystem_row.addWidget(subsystem_alt)

                # Add the row to the target layout
                target_layout.addLayout(subsystem_row)

            # Get the current number of subsystems to maintain indexing consistency
            current_index = len(self.main_filenames)

            # Trigger the file selection and row creation process
            select_multiple_files(section, current_index)

        line = QFrame()
        line.setFrameShape(QFrame.Shape.HLine)
        line.setFrameShadow(QFrame.Shadow.Sunken)
        line.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        layout.addWidget(line)

        # Add relabel mapping section right after the heading
        self.add_relabel_mapping_section(layout, "Standard")

        # Create header labels for main and alternative sections
        main_label = QLabel("Main ORCA Output Files")
        alt_label = QLabel("Alternative ORCA Output Files (Optional)")

        # Styling the labels: equal size, bold text (font 2px smaller than select button text), and height matching select boxes
        main_label.setStyleSheet("background-color: #ebedef; font-size: 11px; padding: 5px; text-align: center; font-weight: bold;")
        alt_label.setStyleSheet("background-color: #ebedef; font-size: 11px; padding: 5px; text-align: center; font-weight: bold;")

        # Set equal size policy for dynamic resizing of both headers
        main_label.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        alt_label.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)

        # Match the height of the select boxes and headings (assuming 30px for the select buttons)
        select_button_height = 30
        main_label.setFixedHeight(select_button_height)
        alt_label.setFixedHeight(select_button_height)

        # Layout for the headers
        header_layout = QHBoxLayout()
        header_layout.addWidget(main_label)  # Main label dynamically resizes
        header_layout.addSpacerItem(header_spacer)  # Fixed 3px spacer between labels
        header_layout.addWidget(alt_label)  # Alt label dynamically resizes

        # Add the header layout to the main layout (below the info button)
        layout.addLayout(header_layout)

        # Supersystem Section (Fixed-size select buttons, expandable populate fields)
        supersystem_layout = QHBoxLayout()

        self.supersystem_main = QLineEdit()
        self.supersystem_alt = QLineEdit()

        # Set placeholders
        self.supersystem_main.setPlaceholderText("Supersystem Main File")
        self.supersystem_alt.setPlaceholderText("Supersystem Alternative File")

        # Set size policy and height to match the heading size
        self.supersystem_main.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.supersystem_alt.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.supersystem_main.setFixedHeight(select_button_height)
        self.supersystem_alt.setFixedHeight(select_button_height)

        # File selection buttons for supersystem files, matching the height of the select and header boxes
        select_supersystem_button = QPushButton("Select Supersystem File")
        select_supersystem_button.setFixedWidth(150)
        select_supersystem_button.setFixedHeight(select_button_height)  # Match height to the header

        select_supersystem_alt_button = QPushButton("Select Supersystem File")
        select_supersystem_alt_button.setFixedWidth(150)
        select_supersystem_alt_button.setFixedHeight(select_button_height)  # Match height to the header

        select_supersystem_button.clicked.connect(lambda: select_file(self.supersystem_main, "main", 0, "Standard", is_supersystem=True))
        select_supersystem_alt_button.clicked.connect(lambda: select_file(self.supersystem_alt, "alt", 0, "Standard", is_supersystem=True))

        # Add widgets to the supersystem layout
        supersystem_layout.addWidget(select_supersystem_button)
        supersystem_layout.addWidget(self.supersystem_main)
        supersystem_layout.addWidget(select_supersystem_alt_button)
        supersystem_layout.addWidget(self.supersystem_alt)

        layout.addLayout(supersystem_layout)

        # Handle manual input for supersystem paths
        self.supersystem_main.editingFinished.connect(lambda: self.update_filename(self.main_filenames, self.normalize_path_gui(self.supersystem_main.text().strip()), 0))
        self.supersystem_alt.editingFinished.connect(lambda: self.update_filename(self.alternative_filenames, self.normalize_path_gui(self.supersystem_alt.text().strip()), 0))

        # Subsystem layout for Standard (below supersystem)
        subsystem_layout = QVBoxLayout()
        layout.addLayout(subsystem_layout)

        # Button to add subsystem files dynamically
        add_button = QPushButton("Add Subsystem Files")
        add_button.setFixedHeight(select_button_height)
        add_button.clicked.connect(lambda: add_subsystems(subsystem_layout, "Standard"))
        layout.addWidget(add_button)

        layout.addSpacing(10)

        # Final spacing and confirm button
        layout.addSpacing(10)
        confirm_button = QPushButton("Confirm Files to Proceed")
        confirm_button.setFixedHeight(40)
        confirm_button.setStyleSheet("QPushButton {background-color: #ebedef;} QPushButton:hover {background-color: paleturquoise;}")
        confirm_button.clicked.connect(self.confirm_standard_nbody_to_proceed)
        layout.addWidget(confirm_button)

        layout.addSpacing(10)

        # Status label for the job process (error message can be copied)
        status_label = QLabel('')
        status_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        status_label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        layout.addWidget(status_label)

        return self.tabs.indexOf(self.nbody_tab)


    def confirm_standard_nbody_to_proceed(self):
        """Confirm the selected files for all sections, ensuring consistent subsystem counts."""

        # Force an update of the supersystem filenames to capture any last-minute changes from the user input fields
        self.main_filenames[0] = self.normalize_path_gui(self.supersystem_main.text().strip() or '')
        self.alternative_filenames[0] = self.normalize_path_gui(self.supersystem_alt.text().strip() or '')

        # Filter out empty entries from each section and validate without overwriting the lists if validation fails
        validated_main, validated_alt = self.filter_empty_rows_and_validate_main(self.main_filenames, self.alternative_filenames, "Standard")
        if validated_main is None:  # Return early if validation failed
            return
        self.main_filenames, self.alternative_filenames = validated_main, validated_alt

        # Count the number of subsystems in each section
        num_subsystems = self.count_valid_subsystems(self.main_filenames)

        # Ensure at least two valid subsystems are provided in each section
        if not (num_subsystems >= 2):
            QMessageBox.warning(self, "Invalid Input", "Please provide one supersystem and at least two subsystem main files in each section.")
            return

        # Prepare the lists of files for confirmation (alternative files can be empty)
        message = "If you want to continue with the following ORCA output files,<br>" \
                  "please press \"OK\". The N-body tab will then be locked to proceed."

        files_overview = f"""<br>
            <b>Main Files:</b> {self.main_filenames}<br>
            <b>Alternative Files:</b> {self.alternative_filenames}<br>
            <b>Relabel Mapping:</b> {self.relabel_mappings.get("Standard", "Not provided")}<br>
        """

        # Create a scrollable confirmation dialog with selectable text
        confirmation_dialog = QDialog(self)
        confirmation_dialog.setWindowTitle("Confirm Files")

        # Remove the "?" help button from the dialog window
        confirmation_dialog.setWindowFlags(confirmation_dialog.windowFlags() & ~Qt.WindowType.WindowContextHelpButtonHint)

        layout = QVBoxLayout()

        # Add a QLabel for the unselectable portion (the instructional message)
        label = QLabel(message)
        layout.addWidget(label)

        # Add QTextEdit for selectable text and set it to read-only
        text_area = QTextEdit()
        text_area.setReadOnly(True)
        text_area.setHtml(files_overview)  # Use setHtml to apply rich text (bold)

        # Allow the text area to dynamically resize with the window
        text_area.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        layout.addWidget(text_area)

        # Add OK and Cancel buttons
        button_box = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        button_box.accepted.connect(confirmation_dialog.accept)
        button_box.rejected.connect(confirmation_dialog.reject)
        layout.addWidget(button_box)

        confirmation_dialog.setLayout(layout)

        # Adjust the dialog size to fit the contents better, and let it dynamically resize
        confirmation_dialog.resize(670, 400)

        # Show the dialog and check user response
        if confirmation_dialog.exec() == QDialog.DialogCode.Accepted:
            self.lock_tab(self.nbody_tab)
            self.switch_to_next_dynamic_tab()


    def add_supersystem_and_subsystems_files(self, layout):
        """Add file selection for the supersystem and subsystems."""
        layout.addWidget(QLabel("Select supersystem and subsystem files:"))

        # Supersystem file selection
        supersystem_main = QLineEdit()
        supersystem_main.setPlaceholderText("Supersystem Main File")
        supersystem_alt = QLineEdit()
        supersystem_alt.setPlaceholderText("Supersystem Alternative File")

        layout.addWidget(supersystem_main)
        layout.addWidget(supersystem_alt)

        # Subsystem file selection logic here (similar to what you had before)
        add_button = QPushButton("Add More Subsystems")
        layout.addWidget(add_button)


    def create_twobody_tab(self):
        """Create the Two-body tab with different content based on CPS and CBS selections."""
        # Check Nbody, CPS, and CBS selections
        nbody_selected = self.nbody_checkbox.isChecked()
        cps_selected = self.cps_checkbox.isChecked()
        cbs_selected = self.cbs_checkbox.isChecked()

        # Create the Two-body tab widget
        self.twobody_tab = QWidget()

        # Create a scroll area for the layout
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)  # Allow the scroll area to resize with the window

        # Create a widget that will be placed inside the scroll area
        scroll_widget = QWidget()

        # Create a layout for the scroll widget
        layout = QVBoxLayout()
        scroll_widget.setLayout(layout)

        # Scenario 1: Both CPS and CBS selected
        if cps_selected and cbs_selected:
            self.add_cps_cbs_twobody_layout(layout)

        # Scenario 2: Only CPS selected
        elif cps_selected and not cbs_selected:
            self.add_cps_twobody_layout(layout)

        # Scenario 3: Only CBS selected
        elif not cps_selected and cbs_selected:
            self.add_cbs_twobody_layout(layout)

        # Scenario 4: Neither CPS nor CBS selected
        else:
            self.add_standard_twobody_layout(layout)

        # Set the scroll widget as the content of the scroll area
        scroll_area.setWidget(scroll_widget)

        # Now, add the scroll area to the Two-body tab instead of setting the layout directly
        layout_wrapper = QVBoxLayout(self.twobody_tab)
        layout_wrapper.addWidget(scroll_area)

        # Insert the Two-body tab before Exit and About
        exit_index = self.tabs.count() - 2  # Assuming Exit is the second last tab
        self.tabs.insertTab(exit_index, self.twobody_tab, "Two-body")
        self.created_tabs.append(self.twobody_tab)

        # Move to the Two-body tab after Home and Nbody confirmation
        self.tabs.setCurrentIndex(self.tabs.indexOf(self.twobody_tab))

        # Return the correct tab index
        return self.tabs.indexOf(self.twobody_tab)


    def reset_button_style(self, move_up_button, move_down_button):
        """Reset the move buttons color back to default."""
        move_up_button.setStyleSheet("background-color: #d3e2f4; color: white; border-radius: 3px;")
        move_down_button.setStyleSheet("background-color: #d3e2f4; color: white; border-radius: 3px;")


    def highlight_up_down_buttons(self, move_up_button, move_down_button):
        """Temporarily change the color of the move up/down buttons."""
        move_up_button.setStyleSheet("background-color: lightblue; color: white;")
        move_down_button.setStyleSheet("background-color: lightblue; color: white;")
        QTimer.singleShot(500, lambda: self.reset_button_style(move_up_button, move_down_button))


    def reset_up_down_button_style(self, move_up_button, move_down_button):
        """Reset the move buttons color back to default."""
        move_up_button.setStyleSheet("background-color: #d3e2f4; color: white; border-radius: 3px;")
        move_down_button.setStyleSheet("background-color: #d3e2f4; color: white; border-radius: 3px;")


    def add_file_selection_row_with_up_down(self, layout, label_text, file_attr, multiple_files=False):
        """Add rows for file selection with an option to move files up/down, remove, and dynamically add new rows."""

        subsection_layout = QVBoxLayout()  # A new layout for this subsection

        # Create the section label with larger font and bold
        label = QLabel(label_text)
        label.setStyleSheet("font-size: 12px; background-color: #d3e2f4; padding: 5px; text-align: center;")
        label.setFixedHeight(30)
        subsection_layout.addWidget(label)

        # Initialize the file list attribute if not already set
        if not hasattr(self, file_attr):
            setattr(self, file_attr, [])  # Store files in this attribute

        file_list = getattr(self, file_attr)
        file_rows = []  # List to store all the file rows and their layouts for this subsection


        def sync_file_list():
            """Synchronize the file_list with the current file_rows order and content."""
            updated_file_list = [file_field.text() for file_field, _ in file_rows]
            setattr(self, file_attr, updated_file_list)


        def add_file_row(file_path=None):
            """Create a new file selection row and add it to the subsection layout."""
            row_layout = QHBoxLayout()

            # File field (with existing file path if any)
            file_field = QLineEdit(file_path if file_path else "")
            file_field.setPlaceholderText("ORCA Output Files")
            file_field.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
            file_field.setFixedHeight(28)

            # Allow manual editing of the file field and update file list accordingly
            file_field.editingFinished.connect(lambda: sync_file_list())

            # Move Up/Down buttons
            move_up_button = QPushButton("▲")
            move_up_button.setFixedWidth(30)
            move_up_button.setFixedHeight(28)
            move_up_button.setStyleSheet("background-color: #d3e2f4; color: white; border-radius: 3px;")
            
            move_down_button = QPushButton("▼")
            move_down_button.setFixedWidth(30)
            move_down_button.setFixedHeight(28)
            move_down_button.setStyleSheet("background-color: #d3e2f4; color: white; border-radius: 3px;")

            remove_button = QPushButton("✖")
            remove_button.setFixedWidth(30)
            remove_button.setFixedHeight(28)
            remove_button.setStyleSheet("background-color: #d3e2f4; color: white; border-radius: 3px;")

            # Add the current file field to the list
            file_rows.append((file_field, row_layout))


            def remove_row():
                """Temporarily change the button color before removing the row."""

                # Change button color to indicate action
                remove_button.setStyleSheet("background-color: red; color: white;")

                # Clear the focus from the current input field to avoid selecting another row
                file_field.clearFocus()

                # Temporarily focus on a non-editable widget (e.g., the main window or a label) to prevent auto-selection
                self.clearFocus()  # Assuming 'self' is your main widget, use it to remove focus from any input field.


                # Delay the actual removal to allow the color change to show
                def delayed_row_removal():
                    # Remove the row and corresponding file from the file list and the layout
                    index = file_rows.index((file_field, row_layout))
                    if index is not None:
                        # Remove from the file list and file_rows list
                        file_list.pop(index)
                        file_rows.pop(index)
                        # Remove the actual row layout from the subsection
                        for i in reversed(range(row_layout.count())):
                            widget = row_layout.itemAt(i).widget()
                            if widget:
                                widget.setParent(None)  # Remove widgets from the row
                        subsection_layout.removeItem(row_layout)  # Remove the layout from the subsection
                        rebuild_subsection_layout()  # Redraw the layout after removal
                        sync_file_list()  # Sync file list after row removal

                    # Ensure no widget gets focus after deletion
                    self.setFocus()  # Focus the main widget, or you could focus another non-editable widget.

                # Use a QTimer to delay the actual deletion for a short time (e.g., 150ms)
                QTimer.singleShot(150, delayed_row_removal)


            def move_up():
                """Move the selected row up."""
                index = file_rows.index((file_field, row_layout))
                if index > 0:
                    # Swap the rows and the corresponding file list elements
                    file_rows[index - 1], file_rows[index] = file_rows[index], file_rows[index - 1]
                    self.highlight_up_down_buttons(move_up_button, move_down_button) # Highlight the buttons after the row is moved
                    rebuild_subsection_layout()  # Redraw the layout after movement
                    sync_file_list()  # Sync file list after movement


            def move_down():
                """Move the selected row down."""
                index = file_rows.index((file_field, row_layout))
                if index < len(file_rows) - 1:
                    # Swap the rows and the corresponding file list elements
                    file_rows[index + 1], file_rows[index] = file_rows[index], file_rows[index + 1]
                    self.highlight_up_down_buttons(move_up_button, move_down_button) # Highlight the buttons after the row is moved
                    rebuild_subsection_layout()  # Redraw the layout after movement
                    sync_file_list()  # Sync file list after movement
        
            # Connect move and remove buttons to their functions
            move_up_button.clicked.connect(move_up)
            move_down_button.clicked.connect(move_down)
            remove_button.clicked.connect(remove_row)

            # Add widgets to the row
            row_layout.addWidget(file_field)
            row_layout.addWidget(move_up_button)
            row_layout.addWidget(move_down_button)
            row_layout.addWidget(remove_button)

            # Add the row to the subsection layout
            subsection_layout.insertLayout(subsection_layout.count() - 1, row_layout)  # Insert before the button row


        def add_file_selection_button():
            """Add a row with only a file selection button (no placeholder field)."""
            row_layout = QHBoxLayout()

            # File selection button
            file_button = QPushButton("Select Files")
            file_button.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)  # Dynamic button resizing with window
            file_button.setFixedHeight(28)


            # File selection logic
            def select_file():
                selected_files, _ = QFileDialog.getOpenFileNames(self, "Select ORCA Output Files", self.last_selected_dir)
                if selected_files:
                    self.last_selected_dir = selected_files[0].rsplit("/", 1)[0]
                    for selected_file in selected_files:
                        # Add each selected file as a new row and update file list
                        file_list.append(self.normalize_path_gui(selected_file))  # Append to file_list first
                        add_file_row(self.normalize_path_gui(selected_file))  # Then add to UI

                    # Ensure the Select Files button stays at the bottom
                    rebuild_subsection_layout()
                    sync_file_list()  # Sync file list after file selection

            # Connect file selection button to the dialog
            file_button.clicked.connect(select_file)

            # Group the button and the line together
            button_layout = QVBoxLayout()
            button_layout.addWidget(file_button)

            # Adding a horizontal line (QFrame) for visual separation
            line = QFrame()
            line.setFrameShape(QFrame.Shape.HLine)
            line.setFrameShadow(QFrame.Shadow.Sunken)
            button_layout.addWidget(line)

            # Return this as a group
            return button_layout


        def rebuild_subsection_layout():
            """Helper function to clear and rebuild the subsection layout after adding/removing/moving rows."""
            # First, clear the layout but keep the label and the Select Files button intact
            for i in reversed(range(subsection_layout.count())):
                item = subsection_layout.itemAt(i)
                widget = item.widget() if item else None
                if widget:
                    widget.setParent(None)
                else:
                    subsection_layout.removeItem(item)

            # Re-add the label and rows in the correct order
            subsection_layout.addWidget(label)  # Ensure the label remains at the top
            for file_field, row_layout in file_rows:
                subsection_layout.addLayout(row_layout)
            subsection_layout.addLayout(file_button_row)  # Ensure the Select Files button group stays at the bottom

        # Initially, add a file selection button row as a group with the horizontal line
        file_button_row = add_file_selection_button()
        subsection_layout.addLayout(file_button_row)

        # Add the subsection layout to the main layout
        layout.addLayout(subsection_layout)


    def add_file_selection_row_without_up_down(self, layout, label_text, file_attr, multiple_files=False):
        """Add rows for Two-body file selection without move up/down buttons."""

        subsection_layout = QVBoxLayout()  # A new layout for this subsection

        # Create the section label with larger font and bold
        label = QLabel(label_text)
        label.setStyleSheet("font-size: 12px; background-color: #d3e2f4; padding: 5px; text-align: center;")
        label.setFixedHeight(30)
        subsection_layout.addWidget(label)

        # Initialize the file list attribute if not already set
        if not hasattr(self, file_attr):
            setattr(self, file_attr, [])  # Store files in this attribute

        file_list = getattr(self, file_attr)
        file_rows = []  # List to store all the file rows and their layouts for this subsection


        def sync_file_list():
            """Synchronize the file_list with the current file_rows order and content."""
            updated_file_list = [file_field.text() for file_field, _ in file_rows]
            setattr(self, file_attr, updated_file_list)


        def add_file_row(file_path=None):
            """Create a new file selection row and add it to the subsection layout."""
            row_layout = QHBoxLayout()

            # File field (with existing file path if any)
            file_field = QLineEdit(file_path if file_path else "")
            file_field.setPlaceholderText("ORCA Output Files")
            file_field.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
            file_field.setFixedHeight(28)

            # Allow manual editing of the file field and update file list accordingly
            file_field.editingFinished.connect(lambda: sync_file_list())

            # Remove button (No Up/Down buttons for Two-body sections)
            remove_button = QPushButton("✖")
            remove_button.setFixedWidth(30)
            remove_button.setFixedHeight(28)
            remove_button.setStyleSheet("background-color: #d3e2f4; color: white; border-radius: 3px;")

            # Add the current file field to the list
            file_rows.append((file_field, row_layout))


            def remove_row():
                """Temporarily change the button color before removing the row."""

                # Change button color to indicate action
                remove_button.setStyleSheet("background-color: red; color: white;")

                # Clear the focus from the current input field to avoid selecting another row
                file_field.clearFocus()

                # Temporarily focus on a non-editable widget (e.g., the main window or a label) to prevent auto-selection
                self.clearFocus()  # Assuming 'self' is your main widget, use it to remove focus from any input field.


                # Delay the actual removal to allow the color change to show
                def delayed_row_removal():
                    # Remove the row and corresponding file from the file list and the layout
                    index = file_rows.index((file_field, row_layout))
                    if index is not None:
                        # Remove from the file list and file_rows list
                        file_list.pop(index)
                        file_rows.pop(index)
                        # Remove the actual row layout from the subsection
                        for i in reversed(range(row_layout.count())):
                            widget = row_layout.itemAt(i).widget()
                            if widget:
                                widget.setParent(None)  # Remove widgets from the row
                        subsection_layout.removeItem(row_layout)  # Remove the layout from the subsection
                        rebuild_subsection_layout()  # Redraw the layout after removal
                        sync_file_list()  # Sync file list after row removal

                    # Ensure no widget gets focus after deletion
                    self.setFocus()  # Focus the main widget, or you could focus another non-editable widget.

                # Use a QTimer to delay the actual deletion for a short time (e.g., 150ms)
                QTimer.singleShot(150, delayed_row_removal)

            # Connect remove button to its function
            remove_button.clicked.connect(remove_row)

            # Add widgets to the row (No Up/Down buttons here)
            row_layout.addWidget(file_field)
            row_layout.addWidget(remove_button)

            # Add the row to the subsection layout
            subsection_layout.insertLayout(subsection_layout.count() - 1, row_layout)  # Insert before the button row


        def add_file_selection_button():
            """Add a row with only a file selection button (no placeholder field)."""
            row_layout = QHBoxLayout()

            # File selection button
            file_button = QPushButton("Select Files")
            file_button.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)  # Dynamic button resizing with window
            file_button.setFixedHeight(28)


            # File selection logic
            def select_file():
                selected_files, _ = QFileDialog.getOpenFileNames(self, "Select ORCA Output Files", self.last_selected_dir)
                if selected_files:
                    self.last_selected_dir = selected_files[0].rsplit("/", 1)[0]
                    for selected_file in selected_files:
                        # Add each selected file as a new row and update file list
                        file_list.append(self.normalize_path_gui(selected_file))  # Append to file_list first
                        add_file_row(self.normalize_path_gui(selected_file))  # Then add to UI

                    # Ensure the Select Files button stays at the bottom
                    rebuild_subsection_layout()
                    sync_file_list()  # Sync file list after file selection

            # Connect file selection button to the dialog
            file_button.clicked.connect(select_file)

            # Group the button and the line together
            button_layout = QVBoxLayout()
            button_layout.addWidget(file_button)

            # Adding a horizontal line (QFrame) for visual separation
            line = QFrame()
            line.setFrameShape(QFrame.Shape.HLine)
            line.setFrameShadow(QFrame.Shadow.Sunken)
            button_layout.addWidget(line)

            # Return this as a group
            return button_layout


        def rebuild_subsection_layout():
            """Helper function to clear and rebuild the subsection layout after adding/removing/moving rows."""
            # First, clear the layout but keep the label and the Select Files button intact
            for i in reversed(range(subsection_layout.count())):
                item = subsection_layout.itemAt(i)
                widget = item.widget() if item else None
                if widget:
                    widget.setParent(None)
                else:
                    subsection_layout.removeItem(item)

            # Re-add the label and rows in the correct order
            subsection_layout.addWidget(label)  # Ensure the label remains at the top
            for file_field, row_layout in file_rows:
                subsection_layout.addLayout(row_layout)
            subsection_layout.addLayout(file_button_row)  # Ensure the Select Files button group stays at the bottom

        # Initially, add a file selection button row as a group with the horizontal line
        file_button_row = add_file_selection_button()
        subsection_layout.addLayout(file_button_row)

        # Add the subsection layout to the main layout
        layout.addLayout(subsection_layout)


    def add_cps_cbs_twobody_layout(self, layout):
        """Create a combined CPS and CBS layout for the Two-body tab, with file selection for both One-body and Two-body sections."""

        # Initialize lists for file paths (One-body and Two-body)
        self.onebody_orcaout_files_SB_LPNO = ['']
        self.onebody_orcaout_files_SB_TPNO = ['']
        self.onebody_orcaout_files_LB_LPNO = ['']
        self.onebody_orcaout_files_LB_TPNO = ['']

        self.twobody_orcaout_files_SB_LPNO = ['']
        self.twobody_orcaout_files_SB_TPNO = ['']
        self.twobody_orcaout_files_LB_LPNO = ['']
        self.twobody_orcaout_files_LB_TPNO = ['']

        # Info button for explanation about ORCA Output Files
        info_button = QPushButton('i')
        info_button.setFixedSize(25, 25)
        info_button.setStyleSheet("border-radius: 12px; background-color: lightgray; color: black; font-weight: bold; font-size: 13px;")


        def show_info():
            QMessageBox.information(self, 'ORCA Output Files Info', 'Please select ORCA output files in its corresponding section for isolated fragments/monomers (One-body) and dimers (Two-body) as well as different computational settings at the supersystem (adduct) geometry.\n\n'
                'Fragments in each dimer must belong to different subsystems used for interaction energy calculation.\n\n'
                'LED data is not required in monomer ORCA output files.\n\n'
                'If N-body LED is selected, the fragment labels used in that section will be applied automatically.\n\n'
                'If N-body LED is not selected, fragments will be labeled according to the order in which one-body ORCA output files are provided in this tab. In this case, if BSSE correction is not applied, only one one-body file per fragment is required. However, if BSSE correction is applied, '
                'two one-body files per dimer are required—one for each fragment. The fragment labels used in the LED interaction energy matrices (ranging from 1 to the number of fragments in the supersystem) depend solely on the order in which you provide the one-body files. '
                'For example, if the one-body file for fragment 1 (from the pair 1–3) is listed second, then fragment 1 will receive label 2 in the LED matrices, unless another one-body file for fragment 1 is provided earlier in the list. '
                'Any additional one-body files associated with the same fragment can appear later in any order without affecting the labeling. Use the up/down buttons or manually edit the file paths to adjust the order of the fragments as needed.\n\n'
                'As an example, if you perform the calculations with cc-pvTZ and cc-pVQZ basis sets and with TCutPNO = 1e-6 and TCutPNO=1e-7, then:\n\n'
                'Smaller Basis Set: cc-pVTZ ; Larger Basis Set: cc-pVQZ\n'
                'Looser PNO: TCutPNO = 1e-6 ; Tighter PNO: TCutPNO = 1e-7')
            
        info_button.clicked.connect(show_info)

        # Info button layout at the top
        info_layout = QHBoxLayout()
        info_layout.addStretch()  # Push the info button to the right
        info_layout.addWidget(info_button)
        layout.addLayout(info_layout)

        # Section for One-body files
        one_body_label = QLabel("ONE-BODY ORCA OUTPUT FILES FOR")
        one_body_label.setStyleSheet("font-size: 14px; font-weight: bold;")
        one_body_label.setFixedHeight(30)
        layout.addWidget(one_body_label)

        # Add file selection rows for One-body sections
        self.add_file_selection_row_with_up_down(layout, "Smaller Basis Set and Looser PNO", 'onebody_orcaout_files_SB_LPNO', multiple_files=True)
        self.add_file_selection_row_with_up_down(layout, "Smaller Basis Set and Tighter PNO", 'onebody_orcaout_files_SB_TPNO', multiple_files=True)
        self.add_file_selection_row_with_up_down(layout, "Larger Basis Set and Looser PNO", 'onebody_orcaout_files_LB_LPNO', multiple_files=True)
        self.add_file_selection_row_with_up_down(layout, "Larger Basis Set and Tighter PNO", 'onebody_orcaout_files_LB_TPNO', multiple_files=True)

        layout.addSpacing(40)

        # Section for Two-body files
        two_body_label = QLabel("TWO-BODY ORCA OUTPUT FILES FOR")
        two_body_label.setStyleSheet("font-size: 14px; font-weight: bold;")
        two_body_label.setFixedHeight(30)
        layout.addWidget(two_body_label)

        # Add file selection rows for Two-body sections
        self.add_file_selection_row_without_up_down(layout, "Smaller Basis Set and Looser PNO", 'twobody_orcaout_files_SB_LPNO', multiple_files=True)
        self.add_file_selection_row_without_up_down(layout, "Smaller Basis Set and Tighter PNO", 'twobody_orcaout_files_SB_TPNO', multiple_files=True)
        self.add_file_selection_row_without_up_down(layout, "Larger Basis Set and Looser PNO", 'twobody_orcaout_files_LB_LPNO', multiple_files=True)
        self.add_file_selection_row_without_up_down(layout, "Larger Basis Set and Tighter PNO", 'twobody_orcaout_files_LB_TPNO', multiple_files=True)

        layout.addSpacing(20)

        # Confirm button
        confirm_button = QPushButton("Confirm Files to Proceed")
        confirm_button.setFixedHeight(40)
        confirm_button.setStyleSheet("QPushButton {background-color: #ebedef;} QPushButton:hover {background-color: paleturquoise;}")
        confirm_button.clicked.connect(self.confirm_cps_cbs_twobody_files)
        layout.addWidget(confirm_button)

        layout.addSpacing(10)

        # Status label for the job process (error message can be copied)
        status_label = QLabel('')
        status_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        status_label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        layout.addWidget(status_label)

        return self.tabs.indexOf(self.twobody_tab)


    def confirm_cps_cbs_twobody_files(self):
        """Confirm the selected files for One-body and Two-body CPS-CBS settings."""

        # Initialize temporary directory and the subdirectories
        self.tmp_dir = os.path.join(self.ledaw_out_root, "tmp")
        self.onebody_orcaout_directory_SB_LPNO = os.path.join(self.tmp_dir, 'onebody_orcaout_directory_SB_LPNO')
        self.onebody_orcaout_directory_SB_TPNO = os.path.join(self.tmp_dir, 'onebody_orcaout_directory_SB_TPNO')
        self.onebody_orcaout_directory_LB_LPNO = os.path.join(self.tmp_dir, 'onebody_orcaout_directory_LB_LPNO')
        self.onebody_orcaout_directory_LB_TPNO = os.path.join(self.tmp_dir, 'onebody_orcaout_directory_LB_TPNO')

        self.twobody_orcaout_directory_SB_LPNO = os.path.join(self.tmp_dir, 'twobody_orcaout_directory_SB_LPNO')
        self.twobody_orcaout_directory_SB_TPNO = os.path.join(self.tmp_dir, 'twobody_orcaout_directory_SB_TPNO')
        self.twobody_orcaout_directory_LB_LPNO = os.path.join(self.tmp_dir, 'twobody_orcaout_directory_LB_LPNO')
        self.twobody_orcaout_directory_LB_TPNO = os.path.join(self.tmp_dir, 'twobody_orcaout_directory_LB_TPNO')

        # Ensure directories exist
        if not os.path.exists(self.tmp_dir):
            os.makedirs(self.tmp_dir)
        if not os.path.exists(self.onebody_orcaout_directory_SB_LPNO):
            os.makedirs(self.onebody_orcaout_directory_SB_LPNO)
        if not os.path.exists(self.onebody_orcaout_directory_SB_TPNO):
            os.makedirs(self.onebody_orcaout_directory_SB_TPNO)
        if not os.path.exists(self.onebody_orcaout_directory_LB_LPNO):
            os.makedirs(self.onebody_orcaout_directory_LB_LPNO)
        if not os.path.exists(self.onebody_orcaout_directory_LB_TPNO):
            os.makedirs(self.onebody_orcaout_directory_LB_TPNO)
        if not os.path.exists(self.twobody_orcaout_directory_SB_LPNO):
            os.makedirs(self.twobody_orcaout_directory_SB_LPNO)
        if not os.path.exists(self.twobody_orcaout_directory_SB_TPNO):
            os.makedirs(self.twobody_orcaout_directory_SB_TPNO)
        if not os.path.exists(self.twobody_orcaout_directory_LB_LPNO):
            os.makedirs(self.twobody_orcaout_directory_LB_LPNO)
        if not os.path.exists(self.twobody_orcaout_directory_LB_TPNO):
            os.makedirs(self.twobody_orcaout_directory_LB_TPNO)

        # Update the file paths and directory paths to capture any last-minute changes from the user input fields and normalize them
        self.onebody_orcaout_files_SB_LPNO = [self.normalize_path_gui(file.strip()) for file in self.onebody_orcaout_files_SB_LPNO]
        self.onebody_orcaout_files_SB_TPNO = [self.normalize_path_gui(file.strip()) for file in self.onebody_orcaout_files_SB_TPNO]
        self.onebody_orcaout_files_LB_LPNO = [self.normalize_path_gui(file.strip()) for file in self.onebody_orcaout_files_LB_LPNO]
        self.onebody_orcaout_files_LB_TPNO = [self.normalize_path_gui(file.strip()) for file in self.onebody_orcaout_files_LB_TPNO]

        self.twobody_orcaout_files_SB_LPNO = [self.normalize_path_gui(file.strip()) for file in self.twobody_orcaout_files_SB_LPNO]
        self.twobody_orcaout_files_SB_TPNO = [self.normalize_path_gui(file.strip()) for file in self.twobody_orcaout_files_SB_TPNO]
        self.twobody_orcaout_files_LB_LPNO = [self.normalize_path_gui(file.strip()) for file in self.twobody_orcaout_files_LB_LPNO]
        self.twobody_orcaout_files_LB_TPNO = [self.normalize_path_gui(file.strip()) for file in self.twobody_orcaout_files_LB_TPNO]

        # Validation: Ensure that the user has provided files for One-body and Two-body
        if not all([self.onebody_orcaout_files_SB_LPNO, self.onebody_orcaout_files_SB_TPNO, 
                    self.onebody_orcaout_files_LB_LPNO, self.onebody_orcaout_files_LB_TPNO,
                    self.twobody_orcaout_files_SB_LPNO, self.twobody_orcaout_files_SB_TPNO, 
                    self.twobody_orcaout_files_LB_LPNO, self.twobody_orcaout_files_LB_TPNO]):
            QMessageBox.warning(self, "Incomplete Input", "Please provide all necessary files before proceeding.")
            return  # Stop further execution if validation fails

        # Check if the number of selected files in each subsection (One-body or two-body) is the same
        num_files_SB_LPNO_onebody = len([f for f in self.onebody_orcaout_files_SB_LPNO if f])
        num_files_SB_TPNO_onebody = len([f for f in self.onebody_orcaout_files_SB_TPNO if f])
        num_files_LB_LPNO_onebody = len([f for f in self.onebody_orcaout_files_LB_LPNO if f])
        num_files_LB_TPNO_onebody = len([f for f in self.onebody_orcaout_files_LB_TPNO if f])

        num_files_SB_LPNO_twobody = len([f for f in self.twobody_orcaout_files_SB_LPNO if f])
        num_files_SB_TPNO_twobody = len([f for f in self.twobody_orcaout_files_SB_TPNO if f])
        num_files_LB_LPNO_twobody = len([f for f in self.twobody_orcaout_files_LB_LPNO if f])
        num_files_LB_TPNO_twobody = len([f for f in self.twobody_orcaout_files_LB_TPNO if f])

        if not (num_files_SB_LPNO_onebody == num_files_SB_TPNO_onebody == num_files_LB_LPNO_onebody == num_files_LB_TPNO_onebody):
            QMessageBox.warning(self, "Inconsistent File Count", "The number of selected files is not the same across all one-body sections. Please ensure they have the same number of files.")
            return  # Stop further execution if the validation fails

        if not (num_files_SB_LPNO_twobody == num_files_SB_TPNO_twobody == num_files_LB_LPNO_twobody == num_files_LB_TPNO_twobody):
            QMessageBox.warning(self, "Inconsistent File Count", "The number of selected files is not the same across all two-body sections. Please ensure they have the same number of files.")
            return  # Stop further execution if the validation fails

        # Any one-body list must have more than one file
        if any(num_files <= 1 for num_files in [num_files_SB_LPNO_onebody, num_files_SB_TPNO_onebody, num_files_LB_LPNO_onebody, num_files_LB_TPNO_onebody]):
            QMessageBox.warning(self, "Invalid One-body File Count", "Each one-body section must have more than one file. Please correct this.")
            return  # Stop further execution if validation fails

        # Ensure the number of two-body files is between n-1 and n*(n-1)/2, where n is the number of one-body files
        n = num_files_SB_LPNO_onebody  # the number of one-body files was already validated above to be the same across all sections

        lower_bound = n - 1 
        upper_bound = n * (n - 1) // 2 

        if not ((lower_bound <= num_files_SB_LPNO_twobody <= upper_bound or num_files_SB_LPNO_twobody == n // 2) and
            (lower_bound <= num_files_SB_TPNO_twobody <= upper_bound or num_files_SB_TPNO_twobody == n // 2) and
            (lower_bound <= num_files_LB_LPNO_twobody <= upper_bound or num_files_LB_LPNO_twobody == n // 2) and
            (lower_bound <= num_files_LB_TPNO_twobody <= upper_bound or num_files_LB_TPNO_twobody == n // 2)):
            QMessageBox.warning(self, "Invalid Two-body File Count", 
                                "The number of two-body files in at least one section is not consistent with that allowed by the number of one-body files. Please correct this.")
            return  # Stop further execution if the validation fails

        # Confirmation dialog
        confirmation_dialog = QDialog(self)
        confirmation_dialog.setWindowTitle("Confirm Files")

        # Remove the "?" help button from the dialog window
        confirmation_dialog.setWindowFlags(confirmation_dialog.windowFlags() & ~Qt.WindowType.WindowContextHelpButtonHint)

        layout = QVBoxLayout()

        message = "If you want to continue with the following ORCA output files and directories,<br>" \
                  "please press \"OK\". The Two-body tab will then be locked to proceed."

        # Add a QLabel for the unselectable portion (the instructional message)
        label = QLabel(message)
        layout.addWidget(label)

        # Add QTextEdit for selectable text and set it to read-only
        orcaout_files_overview = f"""
            <br>
            <b>One-body ORCA Output Files For</b><br>
            <b>Smaller Basis Set and Looser PNO:</b> {self.onebody_orcaout_files_SB_LPNO}<br>
            <b>Smaller Basis Set and Tighter PNO:</b> {self.onebody_orcaout_files_SB_TPNO}<br>
            <b>Larger Basis Set and Looser PNO:</b> {self.onebody_orcaout_files_LB_LPNO}<br>
            <b>Larger Basis Set and Tighter PNO:</b> {self.onebody_orcaout_files_LB_TPNO}<br><br>

            <b>Two-body ORCA Output Files For</b><br>
            <b>Smaller Basis Set and Looser PNO:</b> {self.twobody_orcaout_files_SB_LPNO}<br>
            <b>Smaller Basis Set and Tighter PNO:</b> {self.twobody_orcaout_files_SB_TPNO}<br>
            <b>Larger Basis Set and Looser PNO:</b> {self.twobody_orcaout_files_LB_LPNO}<br>
            <b>Larger Basis Set and Tighter PNO:</b> {self.twobody_orcaout_files_LB_TPNO}<br>
        """

        text_area = QTextEdit()
        text_area.setReadOnly(True)
        text_area.setHtml(orcaout_files_overview)  # Use setHtml to apply rich text (bold)

        # Allow the text area to dynamically resize with the window
        text_area.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        layout.addWidget(text_area)

        # Add OK and Cancel buttons
        button_box = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        button_box.accepted.connect(confirmation_dialog.accept)
        button_box.rejected.connect(confirmation_dialog.reject)
        layout.addWidget(button_box)

        confirmation_dialog.setLayout(layout)

        # Adjust the dialog size to fit the contents better, and let it dynamically resize
        confirmation_dialog.resize(670, 400)

        # Show the dialog and check user response
        if confirmation_dialog.exec() == QDialog.DialogCode.Accepted:
            # Create a temporary directory in the output root, removing it first if it exists
            self.tmp_dir = os.path.join(self.ledaw_out_root, "tmp")
            if os.path.exists(self.tmp_dir):
                shutil.rmtree(self.tmp_dir)
            os.makedirs(self.tmp_dir)

            # Define subdirectories for each one-body and two-body category
            self.onebody_orcaout_directory_SB_LPNO = os.path.join(self.tmp_dir, 'onebody_orcaout_directory_SB_LPNO')
            self.onebody_orcaout_directory_SB_TPNO = os.path.join(self.tmp_dir, 'onebody_orcaout_directory_SB_TPNO')
            self.onebody_orcaout_directory_LB_LPNO = os.path.join(self.tmp_dir, 'onebody_orcaout_directory_LB_LPNO')
            self.onebody_orcaout_directory_LB_TPNO = os.path.join(self.tmp_dir, 'onebody_orcaout_directory_LB_TPNO')

            self.twobody_orcaout_directory_SB_LPNO = os.path.join(self.tmp_dir, 'twobody_orcaout_directory_SB_LPNO')
            self.twobody_orcaout_directory_SB_TPNO = os.path.join(self.tmp_dir, 'twobody_orcaout_directory_SB_TPNO')
            self.twobody_orcaout_directory_LB_LPNO = os.path.join(self.tmp_dir, 'twobody_orcaout_directory_LB_LPNO')
            self.twobody_orcaout_directory_LB_TPNO = os.path.join(self.tmp_dir, 'twobody_orcaout_directory_LB_TPNO')

            # Create directories for one-body outputs
            os.makedirs(self.onebody_orcaout_directory_SB_LPNO, exist_ok=True)
            os.makedirs(self.onebody_orcaout_directory_SB_TPNO, exist_ok=True)
            os.makedirs(self.onebody_orcaout_directory_LB_LPNO, exist_ok=True)
            os.makedirs(self.onebody_orcaout_directory_LB_TPNO, exist_ok=True)

            # Create directories for two-body outputs
            os.makedirs(self.twobody_orcaout_directory_SB_LPNO, exist_ok=True)
            os.makedirs(self.twobody_orcaout_directory_SB_TPNO, exist_ok=True)
            os.makedirs(self.twobody_orcaout_directory_LB_LPNO, exist_ok=True)
            os.makedirs(self.twobody_orcaout_directory_LB_TPNO, exist_ok=True)

            # Copy files to the corresponding directories
            file_mappings = {
                self.onebody_orcaout_directory_SB_LPNO: self.onebody_orcaout_files_SB_LPNO,
                self.onebody_orcaout_directory_SB_TPNO: self.onebody_orcaout_files_SB_TPNO,
                self.onebody_orcaout_directory_LB_LPNO: self.onebody_orcaout_files_LB_LPNO,
                self.onebody_orcaout_directory_LB_TPNO: self.onebody_orcaout_files_LB_TPNO,
                self.twobody_orcaout_directory_SB_LPNO: self.twobody_orcaout_files_SB_LPNO,
                self.twobody_orcaout_directory_SB_TPNO: self.twobody_orcaout_files_SB_TPNO,
                self.twobody_orcaout_directory_LB_LPNO: self.twobody_orcaout_files_LB_LPNO,
                self.twobody_orcaout_directory_LB_TPNO: self.twobody_orcaout_files_LB_TPNO
            }

            for target_dir, file_list in file_mappings.items():
                for file in file_list:
                    if os.path.isfile(file):
                        original_filename = os.path.basename(file)
                        base_name, extension = os.path.splitext(original_filename)
                        destination_path = os.path.join(target_dir, original_filename)
                        counter = 1

                        while os.path.exists(destination_path):
                            new_filename = f"{base_name}_{counter}{extension}"
                            destination_path = os.path.join(target_dir, new_filename)
                            counter += 1

                        shutil.copy(file, destination_path)

            self.lock_tab(self.twobody_tab)
            self.switch_to_next_dynamic_tab()

    def add_cps_twobody_layout(self, layout):
        """Create a CPS layout for the Two-body tab, with file selection for both One-body and Two-body sections."""

        # Initialize lists for file paths (One-body and Two-body)
        self.onebody_orcaout_files_LPNO = ['']
        self.onebody_orcaout_files_TPNO = ['']

        self.twobody_orcaout_files_LPNO = ['']
        self.twobody_orcaout_files_TPNO = ['']

        # Info button for explanation about ORCA Output Files
        info_button = QPushButton('i')
        info_button.setFixedSize(25, 25)
        info_button.setStyleSheet("border-radius: 12px; background-color: lightgray; color: black; font-weight: bold; font-size: 13px;")


        def show_info():
            QMessageBox.information(self, 'ORCA Output Files Info', 'Please select ORCA output files for each setting in its corresponding section for isolated fragments/monomers (One-body) and dimers (Two-body) at the supersystem (adduct) geometry.\n\n'
                'The fragments in the dimers must belong to different subsystems used for interaction energy calculation.\n\n'
                'LED data is not required in monomer ORCA output files.\n\n'
                'If N-body LED is selected, the fragment labels used in that section will be applied automatically.\n\n'
                'If N-body LED is not selected, fragments will be labeled according to the order in which one-body ORCA output files are provided in this tab. If BSSE correction is not applied, only one one-body file per fragment is required. However, if BSSE correction is applied, '
                'two one-body files per dimer are required—one for each fragment. The fragment labels used in the LED interaction energy matrices (ranging from 1 to the number of fragments in the supersystem) depend solely on the order in which you provide the one-body files. '
                'For example, if the one-body file for fragment 1 (from the pair 1–3) is listed second, then fragment 1 will receive label 2 in the LED matrices, unless another one-body file for fragment 1 is provided earlier in the list. '
                'Any additional one-body files associated with the same fragment can appear later in any order without affecting the labeling. Use the up/down buttons or manually edit the file paths to adjust the order of the fragments as needed.\n\n'
                'As an example, if you perform the calculations with TCutPNO= 1e-6 and TCutPNO=1e-7, then:\n\n'
                'Looser PNO: TCutPNO = 1e-6 ; Tighter PNO: TCutPNO = 1e-7')

        info_button.clicked.connect(show_info)

        # Info button layout at the top
        info_layout = QHBoxLayout()
        info_layout.addStretch()  # Push the info button to the right
        info_layout.addWidget(info_button)
        layout.addLayout(info_layout)

        # Section for One-body files
        one_body_label = QLabel("ONE-BODY ORCA OUTPUT FILES FOR")
        one_body_label.setStyleSheet("font-size: 14px;")
        one_body_label.setFixedHeight(30)
        layout.addWidget(one_body_label)

        # Add file selection rows for One-body sections
        self.add_file_selection_row_with_up_down(layout, "Looser PNO", 'onebody_orcaout_files_LPNO', multiple_files=True)
        self.add_file_selection_row_with_up_down(layout, "Tighter PNO", 'onebody_orcaout_files_TPNO', multiple_files=True)

        layout.addSpacing(40)

        # Section for Two-body files
        two_body_label = QLabel("TWO-BODY ORCA OUTPUT FILES FOR")
        two_body_label.setStyleSheet("font-size: 14px;")
        two_body_label.setFixedHeight(30)
        layout.addWidget(two_body_label)

        # Add file selection rows for Two-body sections
        self.add_file_selection_row_without_up_down(layout, "Looser PNO", 'twobody_orcaout_files_LPNO', multiple_files=True)
        self.add_file_selection_row_without_up_down(layout, "Tighter PNO", 'twobody_orcaout_files_TPNO', multiple_files=True)

        layout.addSpacing(20)

        # Confirm button
        confirm_button = QPushButton("Confirm Files to Proceed")
        confirm_button.setFixedHeight(40)
        confirm_button.setStyleSheet("QPushButton {background-color: #ebedef;} QPushButton:hover {background-color: paleturquoise;}")
        confirm_button.clicked.connect(self.confirm_cps_twobody_files)
        layout.addWidget(confirm_button)

        layout.addSpacing(10)

        # Status label for the job process (error message can be copied)
        status_label = QLabel('')
        status_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        status_label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        layout.addWidget(status_label)

        return self.tabs.indexOf(self.twobody_tab)


    def confirm_cps_twobody_files(self):
        """Confirm the selected files for One-body and Two-body CPS settings."""

        # Initialize temporary directory and the subdirectories
        self.tmp_dir = os.path.join(self.ledaw_out_root, "tmp")
        self.onebody_orcaout_directory_LPNO = os.path.join(self.tmp_dir, 'onebody_orcaout_directory_LPNO')
        self.onebody_orcaout_directory_TPNO = os.path.join(self.tmp_dir, 'onebody_orcaout_directory_TPNO')
        self.twobody_orcaout_directory_LPNO = os.path.join(self.tmp_dir, 'twobody_orcaout_directory_LPNO')
        self.twobody_orcaout_directory_TPNO = os.path.join(self.tmp_dir, 'twobody_orcaout_directory_TPNO')

        # Ensure directories exist
        if not os.path.exists(self.tmp_dir):
            os.makedirs(self.tmp_dir)
        if not os.path.exists(self.onebody_orcaout_directory_LPNO):
            os.makedirs(self.onebody_orcaout_directory_LPNO)
        if not os.path.exists(self.onebody_orcaout_directory_TPNO):
            os.makedirs(self.onebody_orcaout_directory_TPNO)
        if not os.path.exists(self.twobody_orcaout_directory_LPNO):
            os.makedirs(self.twobody_orcaout_directory_LPNO)
        if not os.path.exists(self.twobody_orcaout_directory_TPNO):
            os.makedirs(self.twobody_orcaout_directory_TPNO)

        # Update the file paths and directory paths to capture any last-minute changes from the user input fields and normalize them
        self.onebody_orcaout_files_LPNO = [self.normalize_path_gui(file.strip()) for file in self.onebody_orcaout_files_LPNO]
        self.onebody_orcaout_files_TPNO = [self.normalize_path_gui(file.strip()) for file in self.onebody_orcaout_files_TPNO]
        self.twobody_orcaout_files_LPNO = [self.normalize_path_gui(file.strip()) for file in self.twobody_orcaout_files_LPNO]
        self.twobody_orcaout_files_TPNO = [self.normalize_path_gui(file.strip()) for file in self.twobody_orcaout_files_TPNO]

        # Validation: Ensure that the user has provided files for One-body and Two-body
        if not all([self.onebody_orcaout_files_LPNO, self.onebody_orcaout_files_TPNO,
                    self.twobody_orcaout_files_LPNO, self.twobody_orcaout_files_TPNO]):
            QMessageBox.warning(self, "Incomplete Input", "Please provide all necessary files before proceeding.")
            return  # Stop further execution if validation fails

        # Check if the number of selected files in each subsection (One-body or two-body) is the same
        num_files_LPNO_onebody = len([f for f in self.onebody_orcaout_files_LPNO if f])
        num_files_TPNO_onebody = len([f for f in self.onebody_orcaout_files_TPNO if f])
        
        num_files_LPNO_twobody = len([f for f in self.twobody_orcaout_files_LPNO if f])
        num_files_TPNO_twobody = len([f for f in self.twobody_orcaout_files_TPNO if f])
        
        if not (num_files_LPNO_onebody == num_files_TPNO_onebody):
            QMessageBox.warning(self, "Inconsistent File Count", "The number of selected files is not the same across the one-body sections. Please ensure they have the same number of files.")
            return  # Stop further execution if the validation fails

        if not (num_files_LPNO_twobody == num_files_TPNO_twobody):
            QMessageBox.warning(self, "Inconsistent File Count", "The number of selected files is not the same across all two-body sections. Please ensure they have the same number of files.")
            return  # Stop further execution if the validation fails

        # Any one-body list must have more than one file
        if any(num_files <= 1 for num_files in [num_files_LPNO_onebody, num_files_TPNO_onebody]):
            QMessageBox.warning(self, "Invalid One-body File Count", "Each one-body section must have more than one file. Please correct this.")
            return  # Stop further execution if validation fails
        
        # Ensure the number of two-body files is between n-1 and n*(n-1)/2, where n is the number of one-body files
        n = num_files_LPNO_onebody  # the number of one-body files was already validated above to be the same across all sections

        lower_bound = n - 1 
        upper_bound = n * (n - 1) // 2 

        if not ((lower_bound <= num_files_LPNO_twobody <= upper_bound or num_files_LPNO_twobody == n // 2) and
                (lower_bound <= num_files_TPNO_twobody <= upper_bound or num_files_TPNO_twobody == n // 2)):
            QMessageBox.warning(self, "Invalid Two-body File Count", 
                                "The number of two-body files in at least one section is not consistent with that allowed by the number of one-body files. Please correct this.")
            return  # Stop further execution if the validation fails

        # Confirmation dialog
        confirmation_dialog = QDialog(self)
        confirmation_dialog.setWindowTitle("Confirm Files")

        # Remove the "?" help button from the dialog window
        confirmation_dialog.setWindowFlags(confirmation_dialog.windowFlags() & ~Qt.WindowType.WindowContextHelpButtonHint)

        layout = QVBoxLayout()

        message = "If you want to continue with the following ORCA output files and directories,<br>" \
                  "please press \"OK\". The Two-body tab will then be locked to proceed."

        # Add a QLabel for the unselectable portion (the instructional message)
        label = QLabel(message)
        layout.addWidget(label)

        # Add QTextEdit for selectable text and set it to read-only
        orcaout_files_overview = f"""
            <br>
            <b>One-body ORCA Output Files For</b><br>
            <b>Looser PNO:</b> {self.onebody_orcaout_files_LPNO}<br>
            <b>Tighter PNO:</b> {self.onebody_orcaout_files_TPNO}<br><br>

            <b>Two-body ORCA Output Files For</b><br>
            <b>Looser PNO:</b> {self.twobody_orcaout_files_LPNO}<br>
            <b>Tighter PNO:</b> {self.twobody_orcaout_files_TPNO}<br>
        """

        text_area = QTextEdit()
        text_area.setReadOnly(True)
        text_area.setHtml(orcaout_files_overview)  # Use setHtml to apply rich text (bold)

        # Allow the text area to dynamically resize with the window
        text_area.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        layout.addWidget(text_area)

        # Add OK and Cancel buttons
        button_box = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        button_box.accepted.connect(confirmation_dialog.accept)
        button_box.rejected.connect(confirmation_dialog.reject)
        layout.addWidget(button_box)

        confirmation_dialog.setLayout(layout)

        # Adjust the dialog size to fit the contents better, and let it dynamically resize
        confirmation_dialog.resize(670, 400)

        # Show the dialog and check user response
        if confirmation_dialog.exec() == QDialog.DialogCode.Accepted:
            # Define file mappings: map directory names to corresponding file lists
            file_mappings = {
                self.onebody_orcaout_directory_LPNO: self.onebody_orcaout_files_LPNO,
                self.onebody_orcaout_directory_TPNO: self.onebody_orcaout_files_TPNO,
                self.twobody_orcaout_directory_LPNO: self.twobody_orcaout_files_LPNO,
                self.twobody_orcaout_directory_TPNO: self.twobody_orcaout_files_TPNO,
            }

            # Copy the files into the corresponding directories
            for target_dir, file_list in file_mappings.items():
                for file in file_list:
                    if os.path.isfile(file):
                        original_filename = os.path.basename(file)
                        base_name, extension = os.path.splitext(original_filename)
                        destination_path = os.path.join(target_dir, original_filename)
                        counter = 1
    
                        while os.path.exists(destination_path):
                            new_filename = f"{base_name}_{counter}{extension}"
                            destination_path = os.path.join(target_dir, new_filename)
                            counter += 1
    
                        shutil.copy(file, destination_path)

            # Lock the tab and switch to the next dynamic tab
            self.lock_tab(self.twobody_tab)
            self.switch_to_next_dynamic_tab()


    def add_cbs_twobody_layout(self, layout):
        """Create a CBS layout for the Two-body tab, with file selection for both One-body and Two-body sections."""

        # Initialize lists for file paths (One-body and Two-body)
        self.onebody_orcaout_files_SB = ['']
        self.onebody_orcaout_files_LB = ['']

        self.twobody_orcaout_files_SB = ['']
        self.twobody_orcaout_files_LB = ['']

        # Info button for explanation about ORCA Output Files
        info_button = QPushButton('i')
        info_button.setFixedSize(25, 25)
        info_button.setStyleSheet("border-radius: 12px; background-color: lightgray; color: black; font-weight: bold; font-size: 13px;")


        def show_info():
            QMessageBox.information(self, 'ORCA Output Files Info', 'Please select ORCA output files for each setting in its corresponding section for isolated fragments/monomers (One-body) and dimers (Two-body) at the supersystem (adduct) geometry.\n\n'
                'The fragments in the dimers must belong to different subsystems used for interaction energy calculation.\n\n'
                'LED data is not required in monomer ORCA output files.\n\n'
                'If N-body LED is selected, the fragment labels used in that section will be applied automatically.\n\n'
                'If N-body LED is not selected, fragments will be labeled according to the order in which one-body ORCA output files are provided in this tab. If BSSE correction is not applied, only one one-body file per fragment is required. However, if BSSE correction is applied, '
                'two one-body files per dimer are required—one for each fragment. The fragment labels used in the LED interaction energy matrices (ranging from 1 to the number of fragments in the supersystem) depend solely on the order in which you provide the one-body files. '
                'For example, if the one-body file for fragment 1 (from the pair 1–3) is listed second, then fragment 1 will receive label 2 in the LED matrices, unless another one-body file for fragment 1 is provided earlier in the list. '
                'Any additional one-body files associated with the same fragment can appear later in any order without affecting the labeling. Use the up/down buttons or manually edit the file paths to adjust the order of the fragments as needed.\n\n'
                'As an example, if you perform the calculations with cc-pvTZ and cc-pVQZ basis sets, then:\n\n'
                'Smaller Basis Set: cc-pVTZ ; Larger Basis Set: cc-pVQZ')

            
        info_button.clicked.connect(show_info)

        # Info button layout at the top
        info_layout = QHBoxLayout()
        info_layout.addStretch()  # Push the info button to the right
        info_layout.addWidget(info_button)
        layout.addLayout(info_layout)

        # Section for One-body files
        one_body_label = QLabel("ONE-BODY ORCA OUTPUT FILES FOR")
        one_body_label.setStyleSheet("font-size: 14px;")
        one_body_label.setFixedHeight(30)
        layout.addWidget(one_body_label)

        # Add file selection rows for One-body sections
        self.add_file_selection_row_with_up_down(layout, "Smaller Basis Set", 'onebody_orcaout_files_SB', multiple_files=True)
        self.add_file_selection_row_with_up_down(layout, "Larger Basis Set", 'onebody_orcaout_files_LB', multiple_files=True)

        layout.addSpacing(40)

        # Section for Two-body files
        two_body_label = QLabel("TWO-BODY ORCA OUTPUT FILES FOR")
        two_body_label.setStyleSheet("font-size: 14px;")
        two_body_label.setFixedHeight(30)
        layout.addWidget(two_body_label)

        # Add file selection rows for Two-body sections
        self.add_file_selection_row_without_up_down(layout, "Smaller Basis Set", 'twobody_orcaout_files_SB', multiple_files=True)
        self.add_file_selection_row_without_up_down(layout, "Larger Basis Set", 'twobody_orcaout_files_LB', multiple_files=True)

        layout.addSpacing(20)

        # Confirm button
        confirm_button = QPushButton("Confirm Files to Proceed")
        confirm_button.setFixedHeight(40)
        confirm_button.setStyleSheet("QPushButton {background-color: #ebedef;} QPushButton:hover {background-color: paleturquoise;}")
        confirm_button.clicked.connect(self.confirm_cbs_twobody_files)
        layout.addWidget(confirm_button)

        layout.addSpacing(10)

        # Status label for the job process (error message can be copied)
        status_label = QLabel('')
        status_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        status_label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        layout.addWidget(status_label)

        return self.tabs.indexOf(self.twobody_tab)


    def confirm_cbs_twobody_files(self):
        """Confirm the selected files for One-body and Two-body CBS settings."""

        # Initialize temporary directory and the subdirectories
        self.tmp_dir = os.path.join(self.ledaw_out_root, "tmp")
        self.onebody_orcaout_directory_SB = os.path.join(self.tmp_dir, 'onebody_orcaout_directory_SB')
        self.onebody_orcaout_directory_LB = os.path.join(self.tmp_dir, 'onebody_orcaout_directory_LB')
        self.twobody_orcaout_directory_SB = os.path.join(self.tmp_dir, 'twobody_orcaout_directory_SB')
        self.twobody_orcaout_directory_LB = os.path.join(self.tmp_dir, 'twobody_orcaout_directory_LB')

        # Ensure directories exist
        if not os.path.exists(self.tmp_dir):
            os.makedirs(self.tmp_dir)
        if not os.path.exists(self.onebody_orcaout_directory_SB):
            os.makedirs(self.onebody_orcaout_directory_SB)
        if not os.path.exists(self.onebody_orcaout_directory_LB):
            os.makedirs(self.onebody_orcaout_directory_LB)
        if not os.path.exists(self.twobody_orcaout_directory_SB):
            os.makedirs(self.twobody_orcaout_directory_SB)
        if not os.path.exists(self.twobody_orcaout_directory_LB):
            os.makedirs(self.twobody_orcaout_directory_LB)

        # Normalize file paths
        self.onebody_orcaout_files_SB = [self.normalize_path_gui(file.strip()) for file in self.onebody_orcaout_files_SB]
        self.onebody_orcaout_files_LB = [self.normalize_path_gui(file.strip()) for file in self.onebody_orcaout_files_LB]
        self.twobody_orcaout_files_SB = [self.normalize_path_gui(file.strip()) for file in self.twobody_orcaout_files_SB]
        self.twobody_orcaout_files_LB = [self.normalize_path_gui(file.strip()) for file in self.twobody_orcaout_files_LB]

        # Validation: Ensure that all necessary files are provided
        if not all([self.onebody_orcaout_files_SB, self.onebody_orcaout_files_LB, self.twobody_orcaout_files_SB, self.twobody_orcaout_files_LB]):
            QMessageBox.warning(self, "Incomplete Input", "Please provide all necessary files before proceeding.")
            return  # Stop further execution if validation fails

        # Validate consistency of file counts
        num_files_SB_onebody = len([f for f in self.onebody_orcaout_files_SB if f])
        num_files_LB_onebody = len([f for f in self.onebody_orcaout_files_LB if f])
        num_files_SB_twobody = len([f for f in self.twobody_orcaout_files_SB if f])
        num_files_LB_twobody = len([f for f in self.twobody_orcaout_files_LB if f])

        if num_files_SB_onebody != num_files_LB_onebody:
            QMessageBox.warning(self, "Inconsistent File Count", "The number of selected files is not the same across the one-body sections. Please ensure they have the same number of files.")
            return  # Stop further execution if validation fails

        if num_files_SB_twobody != num_files_LB_twobody:
            QMessageBox.warning(self, "Inconsistent File Count", "The number of selected files is not the same across all two-body sections. Please ensure they have the same number of files.")
            return  # Stop further execution if validation fails

        if any(num_files <= 1 for num_files in [num_files_SB_onebody, num_files_LB_onebody]):
            QMessageBox.warning(self, "Invalid One-body File Count", "Each one-body section must have more than one file. Please correct this.")
            return  # Stop further execution if validation fails 

        # Ensure valid number of two-body files
        n = num_files_SB_onebody
        lower_bound = n - 1
        upper_bound = n * (n - 1) // 2

        if not ((lower_bound <= num_files_SB_twobody <= upper_bound or num_files_SB_twobody == n // 2) and
            (lower_bound <= num_files_LB_twobody <= upper_bound or num_files_LB_twobody == n // 2)):
            QMessageBox.warning(self, "Invalid Two-body File Count", 
                                "The number of two-body files in at least one section is not consistent with that allowed by the number of one-body files. Please correct this.")
            return  # Stop further execution if the validation fails

        # Confirmation dialog
        confirmation_dialog = QDialog(self)
        confirmation_dialog.setWindowTitle("Confirm Files")

        # Remove the "?" help button from the dialog window
        confirmation_dialog.setWindowFlags(confirmation_dialog.windowFlags() & ~Qt.WindowType.WindowContextHelpButtonHint)

        layout = QVBoxLayout()

        message = "If you want to continue with the following ORCA output files and directories,<br>" \
                  "please press \"OK\". The Two-body tab will then be locked to proceed."

        # Add a QLabel for the unselectable portion (the instructional message)
        label = QLabel(message)
        layout.addWidget(label)

        # Add QTextEdit for selectable text and set it to read-only
        orcaout_files_overview = f"""<br>
            <b>One-body ORCA Output Files For</b><br>
            <b>Smaller Basis Set:</b> {self.onebody_orcaout_files_SB}<br>
            <b>Larger Basis Set:</b> {self.onebody_orcaout_files_LB}<br><br>

            <b>Two-body ORCA Output Files For</b><br>
            <b>Smaller Basis Set:</b> {self.twobody_orcaout_files_SB}<br>
            <b>Larger Basis Set:</b> {self.twobody_orcaout_files_LB}<br>"""

        text_area = QTextEdit()
        text_area.setReadOnly(True)
        text_area.setHtml(orcaout_files_overview)  # Use setHtml to apply rich text (bold)

        # Allow the text area to dynamically resize with the window
        text_area.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        layout.addWidget(text_area)

        # Add OK and Cancel buttons
        button_box = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        button_box.accepted.connect(confirmation_dialog.accept)
        button_box.rejected.connect(confirmation_dialog.reject)
        layout.addWidget(button_box)

        confirmation_dialog.setLayout(layout)

        # Adjust the dialog size to fit the contents better, and let it dynamically resize
        confirmation_dialog.resize(670, 400)

        # Show the dialog and check user response
        if confirmation_dialog.exec() == QDialog.DialogCode.Accepted:
            # Define file mappings: map directory names to corresponding file lists
            file_mappings = {
                self.onebody_orcaout_directory_SB: self.onebody_orcaout_files_SB,
                self.onebody_orcaout_directory_LB: self.onebody_orcaout_files_LB,
                self.twobody_orcaout_directory_SB: self.twobody_orcaout_files_SB,
                self.twobody_orcaout_directory_LB: self.twobody_orcaout_files_LB,
            }

            # Copy the files into the corresponding directories
            for target_dir, file_list in file_mappings.items():
                for file in file_list:
                    if os.path.isfile(file):
                        original_filename = os.path.basename(file)
                        base_name, extension = os.path.splitext(original_filename)
                        destination_path = os.path.join(target_dir, original_filename)
                        counter = 1
    
                        while os.path.exists(destination_path):
                            new_filename = f"{base_name}_{counter}{extension}"
                            destination_path = os.path.join(target_dir, new_filename)
                            counter += 1
    
                        shutil.copy(file, destination_path)

            # Lock the tab and switch to the next dynamic tab
            self.lock_tab(self.twobody_tab)
            self.switch_to_next_dynamic_tab()


    def add_standard_twobody_layout(self, layout):
        """
        Create a standard layout for the Two-body tab, with file selection for both One-body and Two-body sections.
        """

        # Initialize lists for file paths (One-body and Two-body)
        self.onebody_orcaout_files = ['']
        self.twobody_orcaout_files = ['']

        # Info button for explanation about ORCA Output Files
        info_button = QPushButton('i')
        info_button.setFixedSize(25, 25)
        info_button.setStyleSheet("border-radius: 12px; background-color: lightgray; color: black; font-weight: bold; font-size: 13px;")


        def show_info():
            QMessageBox.information(self, 'ORCA Output Files Info', 'Please select ORCA output files for each setting in its corresponding section for isolated fragments/monomers (One-body) and dimers (Two-body) at the supersystem (adduct) geometry.\n\n'
                'The fragments in the dimers must belong to different subsystems used for interaction energy calculation.\n\n'
                'LED data is not required in monomer ORCA output files.\n\n'
                'If N-body LED is selected, the fragment labels used in that section will be applied automatically.\n\n'
                'If N-body LED is not selected, fragments will be labeled according to the order in which one-body ORCA output files are provided in this tab. If BSSE correction is not applied, only one one-body file per fragment is required. However, if BSSE correction is applied, '
                'two one-body files per dimer are required—one for each fragment. The fragment labels used in the LED interaction energy matrices (ranging from 1 to the number of fragments in the supersystem) depend solely on the order in which you provide the one-body files. '
                'For example, if the one-body file for fragment 1 (from the pair 1–3) is listed second, then fragment 1 will receive label 2 in the LED matrices, unless another one-body file for fragment 1 is provided earlier in the list. '
                'Any additional one-body files associated with the same fragment can appear later in any order without affecting the labeling. Use the up/down buttons or manually edit the file paths to adjust the order of the fragments as needed.\n\n')

            
        info_button.clicked.connect(show_info)

        # Info button layout at the top
        info_layout = QHBoxLayout()
        info_layout.addStretch()  # Push the info button to the right
        info_layout.addWidget(info_button)
        layout.addLayout(info_layout)

        # Add file selection rows for One-body section
        self.add_file_selection_row_with_up_down(layout, "<span style='font-size:11px; font-weight:bold;'>ONE-BODY ORCA OUTPUT FILES</span>", 'onebody_orcaout_files', multiple_files=True)

        layout.addSpacing(40)

        # Add file selection rows for Two-body section
        self.add_file_selection_row_without_up_down(layout, "<span style='font-size:11px; font-weight:bold;'>TWO-BODY ORCA OUTPUT FILES</span>", 'twobody_orcaout_files', multiple_files=True)

        layout.addSpacing(20)

        # Confirm button
        confirm_button = QPushButton("Confirm Files to Proceed")
        confirm_button.setFixedHeight(40)
        confirm_button.setStyleSheet("QPushButton {background-color: #ebedef;} QPushButton:hover {background-color: paleturquoise;}")
        confirm_button.clicked.connect(self.confirm_standard_twobody_files)
        layout.addWidget(confirm_button)

        layout.addSpacing(10)

        # Status label for the job process (error message can be copied)
        status_label = QLabel('')
        status_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        status_label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        layout.addWidget(status_label)

        return self.tabs.indexOf(self.twobody_tab)


    def confirm_standard_twobody_files(self):
        """Confirm the selected files for the standard two-body layout."""

        # Initialize the temporary directory and two-body directory at the beginning
        self.tmp_dir = os.path.join(self.ledaw_out_root, "tmp")
        self.twobody_orcaout_directory = os.path.join(self.tmp_dir, 'twobody_orcaout_directory')
        self.onebody_orcaout_directory = os.path.join(self.tmp_dir, 'onebody_orcaout_directory')

        # Ensure tmp_dir and twobody_orcaout_directory are initialized and exist
        if not os.path.exists(self.tmp_dir):
            os.makedirs(self.tmp_dir)
        if not os.path.exists(self.twobody_orcaout_directory):
            os.makedirs(self.twobody_orcaout_directory)
        if not os.path.exists(self.onebody_orcaout_directory):
            os.makedirs(self.onebody_orcaout_directory)

        # Update the file paths and directory paths to capture any last-minute changes from the user input fields and normalize them
        self.onebody_orcaout_files = [self.normalize_path_gui(file.strip()) for file in self.onebody_orcaout_files]
        self.twobody_orcaout_files = [self.normalize_path_gui(file.strip()) for file in self.twobody_orcaout_files]

        # Validation: Ensure that the user has provided files for One-body and Two-body
        if not all([self.onebody_orcaout_files, self.twobody_orcaout_files]):
            QMessageBox.warning(self, "Incomplete Input", "Please provide all necessary files before proceeding.")
            return  # Stop further execution if validation fails

        # Check the number of selected files in the One-body and Two-body sections
        num_files_onebody = len([f for f in self.onebody_orcaout_files if f])
        num_files_twobody = len([f for f in self.twobody_orcaout_files if f])

        # Any one-body list must have more than one file
        if num_files_onebody <= 1:
            QMessageBox.warning(self, "Invalid One-body File Count", "The one-body section must have more than one file. Please correct this.")
            return  # Stop further execution if validation fails

        # Ensure the number of two-body files is between n-1 and n*(n-1)/2, where n is the number of one-body files
        n = num_files_onebody

        lower_bound = n - 1
        upper_bound = n * (n - 1) // 2

        if not (lower_bound <= num_files_twobody <= upper_bound or num_files_twobody == n // 2):
            QMessageBox.warning(self, "Invalid Two-body File Count",
                                "The number of two-body files is not consistent with that allowed by the number of one-body files. Please correct this.")
            return  # Stop further execution if the validation fails

        # Confirmation dialog
        confirmation_dialog = QDialog(self)
        confirmation_dialog.setWindowTitle("Confirm Files")

        # Remove the "?" help button from the dialog window
        confirmation_dialog.setWindowFlags(confirmation_dialog.windowFlags() & ~Qt.WindowType.WindowContextHelpButtonHint)

        layout = QVBoxLayout()

        message = "If you want to continue with the following ORCA output files and directories,<br>" \
                  "please press \"OK\". The Two-body tab will then be locked to proceed."

        # Add a QLabel for the unselectable portion (the instructional message)
        label = QLabel(message)
        layout.addWidget(label)

        # Add QTextEdit for selectable text and set it to read-only
        orcaout_files_overview = f"""
            <br>
            <b>One-body ORCA Output Files:</b> {self.onebody_orcaout_files}<br><br>

            <b>Two-body ORCA Output Files:</b> {self.twobody_orcaout_files}<br>
        """

        text_area = QTextEdit()
        text_area.setReadOnly(True)
        text_area.setHtml(orcaout_files_overview)  # Use setHtml to apply rich text (bold)

        # Allow the text area to dynamically resize with the window
        text_area.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        layout.addWidget(text_area)

        # Add OK and Cancel buttons
        button_box = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        button_box.accepted.connect(confirmation_dialog.accept)
        button_box.rejected.connect(confirmation_dialog.reject)
        layout.addWidget(button_box)

        confirmation_dialog.setLayout(layout)

        # Adjust the dialog size to fit the contents better, and let it dynamically resize
        confirmation_dialog.resize(670, 400)

        # Show the dialog and check user response
        if confirmation_dialog.exec() == QDialog.DialogCode.Accepted:
            # Define file mappings: map directory names to corresponding file lists
            file_mappings = {self.onebody_orcaout_directory: self.onebody_orcaout_files,
                             self.twobody_orcaout_directory: self.twobody_orcaout_files,}

            # Copy the files into the corresponding directories
            for target_dir, file_list in file_mappings.items():
                for file in file_list:
                    if os.path.isfile(file):
                        original_filename = os.path.basename(file)
                        base_name, extension = os.path.splitext(original_filename)
                        destination_path = os.path.join(target_dir, original_filename)
                        counter = 1
    
                        while os.path.exists(destination_path):
                            new_filename = f"{base_name}_{counter}{extension}"
                            destination_path = os.path.join(target_dir, new_filename)
                            counter += 1
    
                        shutil.copy(file, destination_path)

            # Lock the tab and switch to the next dynamic tab
            self.lock_tab(self.twobody_tab)
            self.switch_to_next_dynamic_tab()


    def create_plot_tab(self):
            """Create the Plot tab with LED heat map parameters."""
            self.plot_tab = QWidget()  # Store the tab in self.plot_tab
            self.plot_tab.setMinimumSize(700, 400)
    
            # Create a scroll area to make the tab content scrollable
            scroll_area = QScrollArea()
            scroll_area.setWidgetResizable(True)
            scroll_area.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)

            # Create a container widget for the content and set its layout
            scroll_content_widget = QWidget()
            scroll_content_widget.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
            layout = QVBoxLayout(scroll_content_widget)
    
            # Create the header for the parameters
            header_layout = QHBoxLayout()
    
            param_header = QLabel("Parameters")
            std_led_header = QLabel("Plot Parameters for Standard LED")
            fp_led_header = QLabel("Plot Parameters for fp-LED")
    
            # Set a consistent height for all header labels
            header_height = 30  # Adjust this value to increase or decrease the height
    
            # Set minimum size for the first header cell (similar to parameter rows)
            min_header_width = 200  # Adjust this value as needed
            param_header.setMinimumWidth(min_header_width)
            std_led_header.setMinimumWidth(min_header_width)
            fp_led_header.setMinimumWidth(min_header_width)
    
            param_header.setMinimumHeight(header_height)
            std_led_header.setMinimumHeight(header_height)
            fp_led_header.setMinimumHeight(header_height)
    
            # Styling the headers
            header_style = "background-color: #ebedef; font-size: 11px; font-weight: bold; padding: 5px;"
            param_header.setStyleSheet(header_style)
            std_led_header.setStyleSheet(header_style)
            fp_led_header.setStyleSheet(header_style)
    
            # Set headers size policy for width flexibility
            param_header.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
            std_led_header.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
            fp_led_header.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
    
            # Add an invisible 'i' button placeholder (taking up space) for the header row
            header_spacer = QSpacerItem(30, header_height, QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)  # Placeholder for alignment
    
            # Add the headers and placeholder
            header_layout.addWidget(param_header)
            header_layout.addWidget(std_led_header)
            header_layout.addWidget(fp_led_header)
            header_layout.addSpacerItem(header_spacer)  # Invisible spacer for alignment
    
            layout.addLayout(header_layout)
    
    
            def create_row(param_text, std_default, fp_default, info_text=None, fig_format_info=False):
                row_layout = QHBoxLayout()
    
                param_input = QLineEdit(param_text)
                param_input.setReadOnly(True)  # Parameter field is unmodifiable
                param_input.setStyleSheet("background-color: lightgrey;")
    
                # Set minimum size for the parameter input (first column)
                param_input.setMinimumWidth(min_header_width)  # Adjust this value based on the text size
                param_input.setSizePolicy(QSizePolicy.Policy.MinimumExpanding, QSizePolicy.Policy.Fixed)  # Allows it to expand horizontally if necessary
    
                # Set consistent size for all input fields
                std_input = QLineEdit(std_default)
                fp_input = QLineEdit(fp_default)
    
                # Set minimum width for the inputs to ensure column alignment
                std_input.setMinimumWidth(min_header_width)
                fp_input.setMinimumWidth(min_header_width)
    
                # Use QSizePolicy.Policy.Expanding for dynamic resizing of the input fields
                std_input.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
                fp_input.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
    
                row_layout.addWidget(param_input)
                row_layout.addWidget(std_input)
                row_layout.addWidget(fp_input)
    
                if info_text or fig_format_info:
                    # Create the visible 'i' button
                    info_button = QPushButton('i')
                    info_button.setFixedSize(24, 24)
                    info_button.setStyleSheet("border-radius: 12px; background-color: lightgray; color: black; font-weight: bold; font-size: 13px;")
                    info_button.clicked.connect(lambda: self.show_info(info_text, fig_format_info=fig_format_info))  # Show the relevant info
                    row_layout.addWidget(info_button)
                else:
                    # If no info text is provided, add an invisible spacer item 
                    spacer = QSpacerItem(30, 30, QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
                    row_layout.addSpacerItem(spacer)
    
                return row_layout, std_input, fp_input
    
            # Add rows for each parameter with 'i' buttons where required
            figure_size_row, self.std_figsize_input, self.fp_figsize_input = create_row(
                "Figure Size", "6, 6", "6, 6",
                info_text='Figure size should be provided as "width, height". Ensure that both width and height are specified.\n\nFor n fragments, a figure size of "n, n" is mostly optimal.'
            )
    
            min_colormap_row, self.std_vmin_input, self.fp_vmin_input = create_row("Minimum of Colormap", "-30", "-30")
            max_colormap_row, self.std_vmax_input, self.fp_vmax_input = create_row("Maximum of Colormap", "30", "30")
            dpi_row, self.std_dpi_input, self.fp_dpi_input = create_row("DPI", "400", "400")
            cutoff_row, self.std_cutoff_input, self.fp_cutoff_input = create_row(
                "Cutoff for Annotation", "None", "None", 
                info_text="If None is selected, all values will be annotated on the heatmap. If a value X is specified, values between -X and X will not be annotated on the heat maps."
            )
            submatrix_row, self.std_submatrix_input, self.fp_submatrix_input = create_row(
                "Enclose a Submatrix with a Black Box", "None", "None", 
                info_text= "To enclose the border of a submatrix with a black box, enter the 0-based indices of the rows and columns "
    					"within this submatrix in the following format: \n\n"
    					"(starting row, starting column), (ending row, ending column) \n\n"
    					"For example, the top-right 3×3 corner of a 6×6 matrix includes rows 0, 1, and 2 (starting row = 0, ending row = 2) "
    					"and columns 3, 4, and 5 (starting column = 3, ending column = 5). "
    					"Thus, in this field you must enter \n\n"
    					"(0, 3), (2, 5) \n\n"
    					"to enclose the top-right 3×3 matrix with a black box. \n\n"
    					"Default: None"
            )
    
            # 'i' button for Figure Format with available options from matplotlib
            fig_format_row, self.std_fig_format_input, self.fp_fig_format_input = create_row(
                "Figure Format", "tif", "tif",
                info_text="Figure format supported by matplotlib: for example, png, jpg, tif, pdf, etc.",
                fig_format_info=True  # Pass this argument to dynamically fetch supported formats
            )
    
            # Add rows to layout
            layout.addLayout(figure_size_row)
            layout.addLayout(min_colormap_row)
            layout.addLayout(max_colormap_row)
            layout.addLayout(fig_format_row)
            layout.addLayout(dpi_row)
            layout.addLayout(cutoff_row)
            layout.addLayout(submatrix_row)
    
            # Initialize the checkbox here to ensure it's an attribute of the class
            self.show_diag_checkbox = QCheckBox("Do you want to demonstrate empty diagonal cells on the fp-LED heat maps?")
            self.show_diag_checkbox.stateChanged.connect(self.toggle_diag_cells_for_fp_led)
            layout.addWidget(self.show_diag_checkbox)
    
            # Similarly for the delete checkbox
            self.delete_old_plots_checkbox = QCheckBox("Do you want to delete heat map directories if they exist from earlier runs?")
            self.delete_old_plots_checkbox.stateChanged.connect(self.toggle_delete_old_plots)
            layout.addWidget(self.delete_old_plots_checkbox)
    
            # Confirm Parameters to Proceed Button
            validate_button = QPushButton("Confirm Parameters to Proceed")
            validate_button.setFixedHeight(40)
            validate_button.setStyleSheet("QPushButton {background-color: #ebedef;} QPushButton:hover {background-color: paleturquoise;}")
            validate_button.clicked.connect(self.confirm_plot_button_clicked)  # Connect to validation and confirmation
            layout.addWidget(validate_button)

            # Set the content widget for the scroll area
            scroll_area.setWidget(scroll_content_widget)
    
            # Create a main layout for the tab and add the scroll area to it
            main_layout = QVBoxLayout(self.plot_tab)
            main_layout.addWidget(scroll_area)
            main_layout.addStretch()
    
            # Insert Plot tab before Exit and About
            exit_index = self.tabs.count() - 2
            self.tabs.insertTab(exit_index, self.plot_tab, "Plot")
            self.created_tabs.append(self.plot_tab)
            return self.tabs.indexOf(self.plot_tab)


    def toggle_diag_cells_for_fp_led(self, checked: bool):
        """Toggle the flag for showing diagonal cells in fp-LED heat maps."""
        self.show_diag_cells_for_fp_led = checked


    def toggle_delete_old_plots(self, checked: bool):
        """Toggle the flag for deleting old heat map directories."""
        self.delete_old_plots = checked


    def show_info(self, text=None, fig_format_info=False):
        """Show information in a message box. Optionally provide figure format info."""
        msg = QMessageBox()
        msg.setWindowTitle("Info")

        try:
            # If running under PyInstaller, use sys._MEIPASS for the path
            if getattr(sys, 'frozen', False):  # Running as a PyInstaller bundle
                base_path = sys._MEIPASS
            else:
                base_path = os.path.dirname(os.path.abspath(sys.argv[0]))

            # Set the icon based on the operating system
            if sys.platform.startswith('win'):
                icon_path = os.path.join(base_path, 'ledaw_package', 'favicon.ico')
            elif sys.platform == 'darwin':  # macOS
                icon_path = os.path.join(base_path, 'ledaw_package', 'favicon.icns')
            else:  # Assume Linux or other OS
                icon_path = os.path.join(base_path, 'ledaw_package', 'favicon.png')

            # Set the message box icon
            msg.setWindowIcon(QIcon(icon_path))

        except Exception as e:
            print(f"Could not set the message box icon: {e}")
            # Continue without crashing if the icon fails to load

        if fig_format_info:
            supported_formats = plt.gcf().canvas.get_supported_filetypes().keys()
            formatted_supported_formats = ', '.join(supported_formats)
            text = f"Supported figure formats for saving plots include: {formatted_supported_formats}"

        msg.setText(text)
        msg.exec()


    def confirm_plot_button_clicked(self):
        """
        Validate and confirm the parameters for the plot.
        Ensure numerical types for figsize, vmin, vmax, dpi, cutoff, and submatrix.
        """
        try:
            # Validation helpers
            def validate_figure_format(fig_format):
                """Check if the figure format is valid based on matplotlib/seaborn supported formats."""
                supported_formats = plt.gcf().canvas.get_supported_filetypes().keys()
                if fig_format.lower() not in supported_formats:
                    raise ValueError(f"Invalid figure format '{fig_format}'. Supported formats: {', '.join(supported_formats)}")


            def validate_cutoff(cutoff_text):
                """Ensure that the cutoff value is either None or a valid float."""
                if cutoff_text.lower() == "none":
                    return None
                try:
                    return float(cutoff_text)
                except ValueError:
                    raise ValueError(f"Invalid cutoff value: {cutoff_text}. Must be 'None' or a valid number.")


            def validate_submatrix(submatrix_text):
                """Ensure the submatrix is either None or a tuple of two tuples."""
                if submatrix_text.lower() == "none":
                    return None
                try:
                    submatrix = eval(submatrix_text)
                    if not (isinstance(submatrix, tuple) and len(submatrix) == 2 and
                            all(isinstance(coords, tuple) and len(coords) == 2 for coords in submatrix)):
                        raise ValueError
                    return submatrix
                except Exception:
                    raise ValueError(f"Invalid submatrix format: {submatrix_text}. Must be None or a tuple of two tuples.")


            def validate_figsize(figsize_text):
                """Ensure that the figure size is a tuple of two elements."""
                try:
                    figsize = tuple(map(float, figsize_text.split(",")))
                    if len(figsize) != 2:
                        raise ValueError
                    return figsize
                except Exception:
                    raise ValueError(f'Invalid figure size format: "{figsize_text}".\n\nIt must consist of two elements (width, height): e.g., (6, 6).')

            # Validate Standard LED parameters
            figsize_std = validate_figsize(self.std_figsize_input.text())
            vmin_std = float(self.std_vmin_input.text())
            vmax_std = float(self.std_vmax_input.text())
            dpi_std = int(self.std_dpi_input.text())
            cutoff_std = validate_cutoff(self.std_cutoff_input.text())
            submatrix_std = validate_submatrix(self.std_submatrix_input.text())
            validate_figure_format(self.std_fig_format_input.text())

            # Validate fp-LED parameters
            figsize_fp = validate_figsize(self.fp_figsize_input.text())
            vmin_fp = float(self.fp_vmin_input.text())
            vmax_fp = float(self.fp_vmax_input.text())
            dpi_fp = int(self.fp_dpi_input.text())
            cutoff_fp = validate_cutoff(self.fp_cutoff_input.text())
            submatrix_fp = validate_submatrix(self.fp_submatrix_input.text())
            validate_figure_format(self.fp_fig_format_input.text())

            # If validation passes, ask for confirmation and lock the tab
            reply = QMessageBox.question(self, "Confirm Parameters",
                                         "Do you want to confirm these settings and lock the Plot tab to proceed?",
                                         QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                                         QMessageBox.StandardButton.No)

            if reply == QMessageBox.StandardButton.Yes:
                # Store validated parameters for Standard LED
                self.plot_params_for_std_led_matrices = {
                    "figsize": figsize_std,
                    "vmin": vmin_std,
                    "vmax": vmax_std,
                    "fig_format": self.std_fig_format_input.text(),
                    "set_dpi": dpi_std,
                    "cutoff_annot": cutoff_std,
                    "submatrix_coords_to_be_highlighted": submatrix_std,
                    "display_heatmap": False
                }

                # Store validated parameters for fp-LED
                self.plot_params_for_fp_led_matrices = {
                    "figsize": figsize_fp,
                    "vmin": vmin_fp,
                    "vmax": vmax_fp,
                    "fig_format": self.fp_fig_format_input.text(),
                    "set_dpi": dpi_fp,
                    "cutoff_annot": cutoff_fp,
                    "submatrix_coords_to_be_highlighted": submatrix_fp,
                    "display_heatmap": False,
                    "show_diag_cells": self.show_diag_cells_for_fp_led
                }

                # Lock the Plot tab after confirming the parameters
                self.lock_tab(self.plot_tab)

                # Automatically switch to the "Run" tab after confirming new parameters
                self.tabs.setCurrentWidget(self.run_tab)
                
                if self.plot_locked:  # If LEDAW has already run
                    self.rerun_plot_button.click()  # Automatically rerun the plot
                    # Re-enable Rerun Plot button for future auto-rerun if needed
                    self.rerun_plot_button.blockSignals(False)  # Allow auto-triggering but keep button unclickable
            
        except ValueError as e:
            # Show validation error in a pop-up dialog
            QMessageBox.warning(self, "Validation Error", f"Validation Error: {str(e)}")


    def create_about_tab(self):
        """Create the About tab."""
        about_tab = QWidget()
        layout = QVBoxLayout()

        # Align the content at the top
        layout.setAlignment(Qt.AlignmentFlag.AlignTop)

        # Create the label with selectable text
        label = QLabel("\nWritten by Ahmet Altun ©\n\nIt is free for academic use. Contact with the author for commercial use.\nFor the full license information, see the LEDAW repository or its manual.\n\nIf you use any part of this code, in addition to original LED, CPS, and CBS studies, please cite:\n1) https://pubs.acs.org/doi/full/10.1021/acs.jcim.5c01561 (J. Chem. Inf. Model. 65, 2025, 8917−8923)\n2) https://github.com/ahmetaltunfatih/LEDAW\n3) https://doi.org/10.1002/anie.202421922 (Angew. Chem. Int. Ed. 64, 2025, e202421922)")

        # Enable text interaction for selecting the text
        label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)

        layout.addWidget(label)

        about_tab.setLayout(layout)
        self.tabs.addTab(about_tab, "About")


    def tab_exists(self, tab_name):
        """Check if a tab with the given name already exists."""
        return any(self.tabs.tabText(i) == tab_name for i in range(self.tabs.count()))


    def create_tabs(self):
        """Create the required tabs based on user selection."""
        created_tabs_indices = []  # Track the indices of created tabs

        # Remove existing N-body, Two-body, and Plot tabs
        self.remove_existing_tabs()

        if self.perform_nbody:
            created_tabs_indices.append(self.create_nbody_tab())
        if self.perform_twobody:
            created_tabs_indices.append(self.create_twobody_tab())
        if self.perform_plot:
            created_tabs_indices.append(self.create_plot_tab())

        return created_tabs_indices


    def add_tabs_in_order(self):
        """Create and switch to the first tab only. Additional tabs will be created one by one after each is locked."""
        # Create the first tab (N-body, Two-body, or Plot) based on user selection.
        if self.perform_nbody:
            self.create_nbody_tab()  # Create N-body tab first
            self.tabs.setCurrentIndex(self.tabs.indexOf(self.nbody_tab))  # Switch to N-body
        elif self.perform_twobody:
            self.create_twobody_tab()  # If N-body is not selected, create Two-body tab first
            self.tabs.setCurrentIndex(self.tabs.indexOf(self.twobody_tab))  # Switch to Two-body
        elif self.perform_plot:
            self.create_plot_tab()  # If neither N-body nor Two-body is selected, create Plot tab first
            self.tabs.setCurrentIndex(self.tabs.indexOf(self.plot_tab))  # Switch to Plot tab directly

        # Disable the Home tab once the first dynamic tab is created
        self.disable_home_tab()

        # Ensure 'Exit' and 'About' tabs are always the last two tabs
        self.ensure_exit_and_about_last()


    def switch_to_next_dynamic_tab(self):
        """Switch to the next dynamically created tab if it exists."""
        current_tab_index = self.tabs.currentIndex()

        # Find the index of the current tab in the created tabs list
        for i, tab in enumerate(self.created_tabs):
            tab_index = self.tabs.indexOf(tab)
            if tab_index == current_tab_index:
                # If this is the current tab, move to the next one
                if i + 1 < len(self.created_tabs):
                    self.tabs.setCurrentIndex(self.tabs.indexOf(self.created_tabs[i + 1]))
                break


    def lock_tab(self, tab_widget):
        """Lock the tab and create the next tab if applicable."""
        if isinstance(tab_widget, QWidget):
            for widget in tab_widget.findChildren(QWidget):
                widget.setEnabled(False)  # Lock the current tab

        # Mark N-body tab as locked, and proceed to create the next tab if requested
        if tab_widget == self.nbody_tab:
            self.nbody_locked = True
            if self.perform_twobody:
                self.create_twobody_tab()  # Create Two-body tab after N-body is locked
            elif self.perform_plot:
                self.create_plot_tab()  # If no Two-body, create Plot tab after N-body is locked

        # Mark Two-body tab as locked, and create Plot tab if requested
        if hasattr(self, 'twobody_tab') and tab_widget == self.twobody_tab:
            self.twobody_locked = True
            if self.perform_plot:
                self.create_plot_tab()  # Create Plot tab after Two-body is locked

        # Mark Plot tab as locked
        if tab_widget == self.plot_tab:
            self.plot_locked = True

        # Check if all dynamic tabs are locked, then create the 'Run' tab
        if (not self.perform_nbody or self.nbody_locked) and \
           (not self.perform_twobody or self.twobody_locked) and \
           (not self.perform_plot or self.plot_locked):
            self.create_run_tab()


    def lock_tabs_in_order(self):
        """Lock all selected tabs in sequence based on user flow."""
        if self.perform_nbody and hasattr(self, 'nbody_tab'):
            self.lock_tab(self.nbody_tab)  # Lock N-body tab first

        if self.perform_twobody and hasattr(self, 'twobody_tab'):
            self.lock_tab(self.twobody_tab)  # Lock Two-body tab if requested

        if self.perform_plot and hasattr(self, 'plot_tab'):
            self.lock_tab(self.plot_tab)  # Lock Plot tab if requested


    def disable_home_tab(self):
        """Disable the interactive elements in the Home tab to prevent modification."""
        for widget in self.home_tab.findChildren(QWidget):
            widget.setEnabled(False)  # Disable all widgets in the Home tab


    def reset_home_tab(self):
        """Re-enable the interactive elements in the Home tab when 'Start Over' is selected."""
        for widget in self.home_tab.findChildren(QWidget):
            widget.setEnabled(True)


    def ensure_exit_and_about_last(self):
        """Ensure 'Exit' and 'About' tabs are always the last two, without recreating them if already present."""
        exit_index = self.tabs.indexOf(self.tabs.findChild(QWidget, "Exit"))
        about_index = self.tabs.indexOf(self.tabs.findChild(QWidget, "About"))

        # If Exit and About tabs are already at the last two positions, do nothing
        if exit_index == self.tabs.count() - 2 and about_index == self.tabs.count() - 1:
            return

        # If not, remove them and append again to ensure they are at the end
        if exit_index != -1:
            exit_tab = self.tabs.widget(exit_index)
            self.tabs.removeTab(exit_index)
            self.tabs.addTab(exit_tab, "Exit")

        if about_index != -1:
            about_tab = self.tabs.widget(about_index)
            self.tabs.removeTab(about_index)
            self.tabs.addTab(about_tab, "About")


    def handle_tab_change(self, index):
        """Handle the tab change event for the Exit tab."""
        if self.tabs.tabText(index) == "Exit":
            self.confirm_exit()


    def confirm_exit(self):
        """Display confirmation dialog for exit."""
        reply = QMessageBox.question(self, "Exit", "Are you sure you want to exit?",
                                     QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                                     QMessageBox.StandardButton.No)

        if reply == QMessageBox.StandardButton.Yes:
            # Properly close the application
            self.close()


    def on_rerun_plot_clicked(self):
        """Handles the UI update for deactivating the Rerun Plot button and creating its copy."""

        # Get the parent layout where the rerun_plot_button exists
        parent_layout = self.rerun_plot_button.parentWidget().layout()

        # Check if the copy button already exists (to avoid adding it multiple times)
        if not hasattr(self, 'rerun_plot_copy_button'):
            # Get the index/position of the rerun_plot_button within the layout
            button_index = parent_layout.indexOf(self.rerun_plot_button)

            # Hide the original button but keep it in the layout
            self.rerun_plot_button.setVisible(False)

            # Create and add the inactive copy button in the same layout position
            self.rerun_plot_copy_button = QPushButton("Rerun Plot", self)
            self.rerun_plot_copy_button.setEnabled(False)  # Disable the button
            self.rerun_plot_copy_button.setFixedHeight(35)

            # Insert the disabled copy button at the same index as the original button
            parent_layout.insertWidget(button_index, self.rerun_plot_copy_button)

            # Apply a custom stylesheet to indicate it's disabled
            self.rerun_plot_copy_button.setStyleSheet("background-color: #ebedef; color: gray;")


    def unlock_plot_tab(self):
        """Unlock the Plot tab to allow user interaction."""
        if hasattr(self, 'plot_tab') and isinstance(self.plot_tab, QWidget):
            # Enable all widgets inside the Plot tab
            for widget in self.plot_tab.findChildren(QWidget):
                widget.setEnabled(True)  # Enable all interactive widgets

            # Enable the Plot tab and make it active
            plot_tab_index = self.tabs.indexOf(self.plot_tab)
            if plot_tab_index != -1:
                self.tabs.setTabEnabled(plot_tab_index, True)  # Enable the Plot tab

            print("Plot tab is unlocked.")

            # Increment the counter each time the Plot tab is unlocked
            self.plot_unlock_count += 1

            # Replace Rerun Plot button only after the first unlock
            if self.plot_unlock_count > 1:
                self.on_rerun_plot_clicked()  # Replace the Rerun Plot button with its inactive copy


    def delete_tmp_files(self):
        """Delete the temporary directory created during the LEDAW run."""

        # Ensure self.ledaw_out_root is not None
        if self.ledaw_out_root is None:
            print("LEDAW output root directory is not set. No temporary files to delete.")
            return

        tmp_directory = os.path.join(self.ledaw_out_root, "tmp")

        # Check if the tmp directory exists and remove it
        if os.path.exists(tmp_directory):
            try:
                shutil.rmtree(tmp_directory) # tmp directory has been deleted.
            except Exception as e:
                print(f"Error while deleting temporary directory: {e}")


    def create_run_tab(self):
        """Create the 'Run' tab after the last dynamically created tab is locked."""

        # Check if Run tab already exists to prevent duplicates
        if hasattr(self, 'run_tab') and isinstance(self.run_tab, QWidget):
            return  # If Run tab exists, don't create it again

        self.setGeometry(100, 100, 600, 400)

        run_tab = QWidget()
        layout = QVBoxLayout()

        # Information button ('i')
        info_button = QPushButton('i')
        info_button.setFixedSize(24, 24)
        info_button.setStyleSheet("border-radius: 12px; background-color: lightgray; font-weight: bold; font-size: 13px;")

        # Info text for the button
        info_text = (
            '"Run LEDAW" performs LED analyses based on the inputs in the previous tabs.\n\n'
            'When the LEDAW run is completed, check the plot directories. If the current figure settings '
            'are not optimal, adjust and then confirm them in the "Plot" tab to rerun the plot job.\n\n'
            'A running job can be cancelled at any time by pressing the "Cancel Run" button.'
        )

        info_button.clicked.connect(lambda: self.show_info(info_text))

        # Align the 'i' button to the right
        info_layout = QHBoxLayout()
        info_layout.addStretch()  # Pushes the button to the right
        info_layout.addWidget(info_button)

        # Add info layout to the main layout
        layout.addLayout(info_layout)

        # Run LEDAW, Rerun Plot, Cancel Run buttons
        self.run_ledaw_button = QPushButton("Run LEDAW")
        self.rerun_plot_button = QPushButton("Rerun Plot")
        self.cancel_run_button = QPushButton("Cancel Run")

        # Set larger heights for the buttons
        self.run_ledaw_button.setFixedHeight(35)
        self.rerun_plot_button.setFixedHeight(35)
        self.cancel_run_button.setFixedHeight(35)

        # Set style of the buttons
        self.run_ledaw_button.setStyleSheet("QPushButton {background-color: #ebedef;} QPushButton:hover {background-color: paleturquoise;}")
        self.rerun_plot_button.setStyleSheet("QPushButton {background-color: #ebedef;} QPushButton:hover {background-color: paleturquoise;}")
        self.cancel_run_button.setStyleSheet("QPushButton {background-color: #ebedef;} QPushButton:hover {background-color: paleturquoise;}")

        # Initially, only Run LEDAW is active, Rerun Plot will be activated after LEDAW job finishes if Plot tab exists
        self.run_ledaw_button.setEnabled(True)
        self.rerun_plot_button.setEnabled(False)  # Initially disabled
        self.cancel_run_button.setEnabled(False)

        # Disable the Cancel Run button from receiving focus entirely
        self.cancel_run_button.setFocusPolicy(Qt.FocusPolicy.NoFocus)

        # Create status labels for each button
        run_ledaw_status_label = QLabel("")
        rerun_plot_status_label = QLabel("")
        cancel_run_status_label = QLabel("")

        # Align the status label text to the center of the row
        run_ledaw_status_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        rerun_plot_status_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        cancel_run_status_label.setAlignment(Qt.AlignmentFlag.AlignCenter)

        # Add buttons to the layout along with status labels
        layout.addWidget(self.run_ledaw_button)
        layout.addWidget(run_ledaw_status_label)
        layout.addWidget(self.rerun_plot_button)
        layout.addWidget(rerun_plot_status_label)
        layout.addWidget(self.cancel_run_button)
        layout.addWidget(cancel_run_status_label)

        # Set the layout for the run tab
        run_tab.setLayout(layout)

        # Insert Run tab after the last dynamically created tab
        last_dynamic_tab_index = self.tabs.indexOf(self.created_tabs[-1]) if self.created_tabs else self.tabs.count() - 2
        self.tabs.insertTab(last_dynamic_tab_index + 1, run_tab, "Run")

        # Set 'Run' tab as current tab
        self.tabs.setCurrentWidget(run_tab)

        # Add 'Run' tab to the list of dynamically created tabs for future deletion
        self.created_tabs.append(run_tab)

        # Store the run_tab for future reference
        self.run_tab = run_tab


        # Function to fetch updated parameters from Plot tab
        def fetch_plot_parameters():
            """Fetch the updated plot parameters from the Plot tab."""
            if hasattr(self, 'plot_tab') and self.plot_tab.isEnabled():
                try:
                    # Fetch updated plot parameters (assuming you have input fields in the Plot tab)
                    figsize_std_str = self.std_figsize_input.text()
                    vmin_std = float(self.std_vmin_input.text())
                    vmax_std = float(self.std_vmax_input.text())
                    dpi_std = int(self.std_dpi_input.text())
                    cutoff_std = self.std_cutoff_input.text()
                    submatrix_std = self.std_submatrix_input.text()

                    # Similarly fetch for fp-LED if applicable
                    figsize_fp_str = self.fp_figsize_input.text()
                    vmin_fp = float(self.fp_vmin_input.text())
                    vmax_fp = float(self.fp_vmax_input.text())
                    dpi_fp = int(self.fp_dpi_input.text())
                    cutoff_fp = self.fp_cutoff_input.text()
                    submatrix_fp = self.fp_submatrix_input.text()

                    # Convert string figsize (e.g., "10, 8") into a tuple (10, 8)
                    figsize_std = tuple(map(float, figsize_std_str.split(',')))
                    figsize_fp = tuple(map(float, figsize_fp_str.split(',')))

                    # Store these updated parameters
                    self.plot_params_for_std_led_matrices.update({
                        "figsize": figsize_std,
                        "vmin": vmin_std,
                        "vmax": vmax_std,
                        "dpi": dpi_std,
                        "cutoff": cutoff_std,
                        "submatrix": submatrix_std
                    })

                    self.plot_params_for_fp_led_matrices.update({
                        "figsize": figsize_fp,
                        "vmin": vmin_fp,
                        "vmax": vmax_fp,
                        "dpi": dpi_fp,
                        "cutoff": cutoff_fp,
                        "submatrix": submatrix_fp
                    })

                    print("Plot parameters have been updated successfully.")

                except ValueError as e:
                    print(f"Error while fetching plot parameters: {e}")


        def run_ledaw():
            print("Running LEDAW...")

            # Start with the "In progress..." message
            run_ledaw_status_label.setText("In progress...")

            # Set up blinking effect
            self.blinking = True  # Start the blinking flag
            self.blink_timer = QTimer(self)  # Set up the timer
            self.blink_timer.timeout.connect(toggle_blink)  # Connect the timer to the blinking toggle function
            self.blink_timer.start(500)  # Blink every 0.5 second

            # Disable Run LEDAW button while the job is running, enable Cancel Run button
            self.run_ledaw_button.setEnabled(False)
            self.cancel_run_button.setEnabled(True)

            # Start the worker thread to run the LEDAW job in the background
            self.worker = JobWorker(self)
            self.worker.status_signal.connect(update_ledaw_run_status)  # Connect worker status signal to update function
            self.worker.start()

            # Reset the cancellation flag
            self.job_cancelled = False
            self.running_job = "ledaw"  # Mark the job as running


        def update_ledaw_run_status(message):
            """Update the status label based on the job progress and handle button states."""
            run_ledaw_status_label.setText(message)

            if "Completed" in message:
                self.blink_timer.stop()  # Stop the blinking effect
                self.run_ledaw_button.setEnabled(False)  # Keep Run LEDAW disabled after completion
                self.cancel_run_button.setEnabled(False)  # Disable Cancel Run button

                # Enable Rerun Plot if plotting is required
                if self.perform_plot:
                    self.rerun_plot_button.setEnabled(True)
                    self.rerun_plot_button.setFocus()
                    self.unlock_plot_tab()
                    self.plot_locked = False
                    self.running_job = None 

            elif "Cancelled" in message:
                self.blink_timer.stop()  # Stop the blinking effect if job was cancelled
                self.run_ledaw_button.setEnabled(True)  # Allow running LEDAW again
                self.cancel_run_button.setEnabled(False)

            elif "Failed" in message:
                self.blink_timer.stop()  # Stop blinking on failure
                self.run_ledaw_button.setEnabled(True)  # Allow resubmission in case of failure
                self.cancel_run_button.setEnabled(False)


        def toggle_blink():
            """Toggle between 'In progress...' and an empty string to simulate blinking."""
            current_text = run_ledaw_status_label.text()

            if current_text == "In progress...":
                run_ledaw_status_label.setText("")  # Clear the label to simulate the blink
            else:
                run_ledaw_status_label.setText("In progress...")  # Show "In progress..." again


        def ledaw_finished():
            """Called when the LEDAW run is finished."""

            # Stop blinking effect
            if hasattr(self, 'blink_timer') and self.blink_timer.isActive():
                self.blink_timer.stop()

            # Check if the job was cancelled
            if self.job_cancelled:
                run_ledaw_status_label.setText("LEDAW job cancelled.")
                self.run_ledaw_button.setEnabled(True)  # Allow rerunning LEDAW if job was cancelled
                self.cancel_run_button.setEnabled(False)  # Disable Cancel Run button
                return  # Skip the normal completion flow if the job was cancelled

            # Job completed successfully
            run_ledaw_status_label.setText("LEDAW job completed successfully.")

            # Disable the Run LEDAW button as no more LEDAW job submission is allowed
            self.run_ledaw_button.setEnabled(False)

            # Enable Rerun Plot button if plotting is required
            if self.perform_plot:
                self.rerun_plot_button.setEnabled(True)  # Activate Rerun Plot button
                self.rerun_plot_button.setFocus()  # Focus on the Rerun Plot button


        def rerun_plot():
            
            self.delete_tmp_files()
            
            if not self.plot_locked:
                print("Please first adjust and then confirm plot parameters in the Plot tab")
                rerun_plot_status_label.setText("Please first adjust and then confirm plot parameters in the Plot tab.")
                run_ledaw_status_label.setText("")  # Clear the other labels
                cancel_run_status_label.setText("")  # Clear the other labels

            else:
                print("Rerunning Plot...")

                rerun_plot_status_label.setText("In progress...")  # Set initial text to "In progress..."
                run_ledaw_status_label.setText("")  # Clear the other labels
                cancel_run_status_label.setText("")  # Clear the other labels

                # Fetch updated parameters from Plot tab before rerunning plot
                fetch_plot_parameters()

                self.rerun_plot_button.setFocus()  # Ensure focus stays on this button
                self.rerun_plot_button.setEnabled(False)  # Disable the button during the job
                self.cancel_run_button.setEnabled(True)  # Activate the Cancel Run button when plot is running

                # Set up blinking effect for "In progress..."
                self.blinking = True  # Track blinking state
                self.blink_timer = QTimer(self)
                self.blink_timer.timeout.connect(toggle_blink_plot)  # Connect the timer to the blinking function
                self.blink_timer.start(500)  # Blink every 0.5 second

                # Create the worker thread and connect the status signal to the label update
                self.worker = PlotJobWorker(self)
                self.worker.status_signal.connect(update_rerun_plot_status)  # Update the status label based on worker signals
                self.worker.start()

                # Reset cancellation flag and track running job
                self.job_cancelled = False  
                self.running_job = "plot"  


        def update_rerun_plot_status(message):
            """Handles updating the status and stops blinking when the plot finishes."""
            rerun_plot_status_label.setText(message)

            if "Completed" in message or "Failed" in message:
                self.blink_timer.stop()  # Stop blinking once the plot job is done
                self.blinking = False  # Ensure blinking state is reset
                self.rerun_plot_button.setEnabled(True)  # Re-enable Rerun Plot button
                self.cancel_run_button.setEnabled(False)  # Disable Cancel Run button when plot finishes
                self.unlock_plot_tab()
                self.running_job = None  # Reset running job status
            elif "Cancelled" in message:
                self.blink_timer.stop()  # Stop blinking if the job is canceled
                self.blinking = False
                self.rerun_plot_button.setEnabled(True)
                self.cancel_run_button.setEnabled(False)
                self.unlock_plot_tab()
                self.running_job = None  # Reset running job status


        def toggle_blink_plot():
            """Toggle between 'In progress...' and an empty string to simulate blinking."""
            if self.blinking:
                rerun_plot_status_label.setText("In progress...")
            else:
                rerun_plot_status_label.setText("")
            self.blinking = not self.blinking  # Toggle the blinking state


        def plot_finished():
            """Called when the Plot run is finished."""
            if self.job_cancelled:
                rerun_plot_status_label.setText(
                    "Cancelled\n\nTo rerun the job, arrange the plot parameters and press 'Confirm Parameters to Proceed' in the Plot tab."
                )
                return  # Skip the normal completion flow if the job was canceled

            print("Plot rerun completed.")

            # Replace the Rerun Plot button with its inactive copy
            self.on_rerun_plot_clicked()  # Handles UI update for button state

            self.cancel_run_button.setEnabled(False)  # Disable Cancel Run button
            self.rerun_plot_button.setFocus()  # Ensure focus on Rerun Plot button (if needed)

            # Reset the running job status
            self.running_job = None


        def confirm_plot_tab():
            """Confirm the plot tab and automatically start the plot job."""
            print("Plot tab confirmed, starting plot job...")

            # Simulate confirmation of the plot tab
            self.plot_locked = True  # Lock the plot to prevent further changes
            rerun_plot()  # Automatically run the plot job


        def cancel_run():
            print("Run cancelled.")

            # Stop the blinking effect
            if hasattr(self, 'blink_timer') and self.blink_timer.isActive():
                self.blink_timer.stop()

            # Clear status labels and update the cancel message
            cancel_run_status_label.setText("")
            run_ledaw_status_label.setText("")  # Clear other labels
            rerun_plot_status_label.setText("")  # Clear other labels
            self.job_cancelled = True  # Mark the job as cancelled

            # Handle button states based on the running job
            if self.running_job == "ledaw":
                # If cancelling LEDAW, reactivate the Run LEDAW button and do not unlock the Plot tab
                self.run_ledaw_button.setEnabled(True)  # Reactivate Run LEDAW button if needed
                self.run_ledaw_button.setFocus()  # Ensure focus goes back to the Run LEDAW button
                print("LEDAW job cancelled. Plot tab will remain locked.")
            elif self.running_job == "plot":
                # If cancelling Rerun Plot, unlock the Plot tab and keep Run LEDAW inactive
                self.rerun_plot_button.setEnabled(True)  # Allow Rerun Plot to be re-run
                self.rerun_plot_button.setFocus()  # Focus back on Rerun Plot button
                self.unlock_plot_tab()  # Unlock the Plot tab so the user can adjust parameters
                print("Plot job cancelled. Plot tab unlocked.")

            # Deactivate Cancel Run button
            self.cancel_run_button.setEnabled(False)

            # Cancel the running job in the worker thread
            if self.worker and self.worker.isRunning():
                self.worker.cancel()  # This stops the worker thread
                print(f"{self.running_job.capitalize()} job cancelled.")
            else:
                print("No job was running.")

            # Reset the running job flag
            self.running_job = None

            
        # Connect buttons to their respective functions
        self.run_ledaw_button.clicked.connect(run_ledaw)
        self.rerun_plot_button.clicked.connect(rerun_plot)
        self.cancel_run_button.clicked.connect(cancel_run)

        # Automatically trigger the appropriate button
        if self.run_ledaw_button.isEnabled():
            self.run_ledaw_button.click()  # Automatically trigger Run LEDAW
        elif self.rerun_plot_button.isEnabled():
            self.rerun_plot_button.click()  # Automatically trigger Rerun Plot     
