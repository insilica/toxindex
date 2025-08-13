import Layout from "./components/Layout";
import LoginForm from "./components/login/LoginForm";
import RegisterForm from "./components/login/RegisterForm";
import VerifyPage from "./components/login/VerifyPage";
import ProtectedRoute from "./components/ProtectedRoute";
// import CreateEnvironment from "./components/CreateEnvironment";
import Dashboard from "./components/Dashboard";
import SettingsEnvironments from "./components/settings/SettingsEnvironments";
import TermsPrivacy from "./components/TermsPrivacy";
import Settings from "./components/settings/Settings";
import SettingsDataControls from "./components/settings/SettingsDataControls";
import CreateEnvironmentSettings from "./components/settings/CreateEnvironmentSettings";
import { BrowserRouter as Router, Routes, Route } from "react-router-dom";
import EnvironmentDetails from "./components/EnvironmentDetails";

import ChatSession from "./components/ChatSession";
import TaskDetail from './components/TaskDetail';
import UserProfile from './components/UserProfile';
import ForgotPasswordPage from './components/login/ForgotPasswordPage';
import ResetPasswordPage from './components/login/ResetPasswordPage';
import AdminRoute from './components/admin/AdminRoute';
import UserGroupManager from './components/admin/UserGroupManager';
import { EnvironmentProvider } from "./context/EnvironmentContext";
import { ChatSessionProvider } from "./context/ChatSessionContext";
import { ModelProvider } from "./context/ModelContext";
import { SocketProvider } from "./context/SocketContext";
import { SessionProvider, useSession } from "./context/SessionContext";
import { SessionTimeoutNotification } from "./components/SessionTimeoutNotification";
import { useEffect } from "react";

// Component for protected routes with session management
function ProtectedRoutes() {
  console.log('[App] ProtectedRoutes component rendering');
  
  // Get session state from context
  const {
    extendSession,
    logoutNow,
    showWarning,
    timeRemaining
  } = useSession();

  console.log('[App] ProtectedRoutes session state:', { showWarning, timeRemaining });

  return (
    <>
      <SessionTimeoutNotification
        isVisible={showWarning}
        onExtendSession={extendSession}
        onLogout={logoutNow}
        timeRemaining={timeRemaining}
      />
      <Routes>
        <Route
          path="/"
          element={
            <ProtectedRoute>
              <Layout>
                <Dashboard />
              </Layout>
            </ProtectedRoute>
          }
        />
        <Route
          path="/settings/general"
          element={
            <ProtectedRoute>
              <Layout>
                <Settings/>
              </Layout>
            </ProtectedRoute>
          }
        />
        <Route
          path="/settings/data-controls"
          element={
            <ProtectedRoute>
              <Layout>
                <SettingsDataControls />
              </Layout>
            </ProtectedRoute>
          }
        />
        <Route
          path="/settings/environments"
          element={
            <ProtectedRoute>
              <Layout>
                <SettingsEnvironments />
              </Layout>
            </ProtectedRoute>
          }
        />
        <Route
          path="/settings/environments/create"
          element={
            <ProtectedRoute>
              <Layout>
                <CreateEnvironmentSettings />
              </Layout>
            </ProtectedRoute>
          }
        />
        <Route
          path="/environments/details"
          element={
            <ProtectedRoute>
              <Layout>
                <EnvironmentDetails />
              </Layout>
            </ProtectedRoute>
          }
        />
        <Route
          path="/chat/session/:sessionId"
          element={
            <ProtectedRoute>
              <Layout>
                <ChatSession />
              </Layout>
            </ProtectedRoute>
          }
        />
        <Route path="/task/:task_id" element={<TaskDetail />} />
        <Route path="/user/:user_id" element={<UserProfile />} />

        {/* Admin routes */}
        <Route
          path="/admin/users"
          element={
            <ProtectedRoute>
              <Layout>
                <AdminRoute>
                  <UserGroupManager />
                </AdminRoute>
              </Layout>
            </ProtectedRoute>
          }
        />
      </Routes>
    </>
  );
}

function AppContent() {
  console.log('[App] AppContent component rendering');
  
  return (
    <Routes>
      {/* Public routes: no session management */}
      <Route path="/login" element={<LoginForm />} />
      <Route path="/register" element={<RegisterForm />} />
      <Route path="/verify/:token" element={<VerifyPage />} />
      <Route path="/policies/terms-of-use/" element={<TermsPrivacy />} />
      <Route path="/policies/privacy-policy/" element={<TermsPrivacy />} />
      <Route path="/forgot_password" element={<ForgotPasswordPage />} />
      <Route path="/reset_password/:token" element={<ResetPasswordPage />} />

      {/* Protected routes: with session management */}
      <Route path="/*" element={
        <SessionProvider onSessionExpired={() => {
          console.log('[App] Session expired, redirecting to login');
          // Clear any cached data
          localStorage.clear();
          sessionStorage.clear();
          window.location.href = '/login';
        }}>
          <ProtectedRoutes />
        </SessionProvider>
      } />
    </Routes>
  );
}

function App() {
  console.log('[App] App component mounting');

  useEffect(() => {
    console.log('[App] App component mounted');
    return () => {
      console.log('[App] App component unmounting');
    };
  }, []);

  console.log('[App] App component rendering');

  return (
    <div className="w-screen h-screen min-h-screen min-w-full">
      <SocketProvider>
        <ModelProvider>
          <ChatSessionProvider>
            <EnvironmentProvider>
              <Router>
                <AppContent />
              </Router>
            </EnvironmentProvider>
          </ChatSessionProvider>
        </ModelProvider>
      </SocketProvider>
    </div>
  );
}

export default App;
